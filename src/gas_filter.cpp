#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <RcppParallel.h>
#include <cstring>
#include <vector>
using namespace Rcpp;

// The score-driven recursion and the derivative of its state, propagated
// together. The callbacks into R -- the score and the curvature belong to
// whatever distribution the model carries -- remain the default; where the
// caller supplies a FAST context (piano_parallel.txt, section 2a), the
// family's and the link's scalar C entry points are resolved once per call
// with R_GetCCallable and the loop runs without touching the R API, which
// is the precondition for the threads over groups below, not only a speed.
//
// The fast step MIRRORS distrib_kernel()'s composition expression by
// expression -- theta_k = clamp(linkinv(e)), s = l_k * h', c = l_kk h'^2 +
// l_k h'' -- and a twin test holds the two routes bit-identical. A family
// or a link the registries do not cover makes the context inert and the
// callbacks run as before: coverage is a speed property, never a
// correctness one.
//
// Every index arriving from R is 1-based and is used to address the
// caller's rows, which is why `row` is passed to the callbacks unchanged.

namespace {

typedef int (*lf7_id_t)(const char*);
typedef void (*lf7_inv12_t)(int, double, double*, double*, double*);
typedef double (*lf7_clamp_t)(double, double, double);
typedef int (*d7_id_t)(const char*);
typedef void (*d7_sc_t)(int, int, double, const double*, double*);

constexpr int kMaxPar = 16;
// below this many groups the opening of the parallel region is not worth
// its cost; measured on the panel benchmarks and held to the plan's
// asymmetry argument (a factor of two either way costs almost nothing)
constexpr int kMinGroupsPar = 8;

struct FastCtx {
    bool ok = false;
    int fam = -1, link = -1, k = 0, npar = 0;
    double lwr = 0.0, upr = 0.0;
    const double* y = nullptr;
    const double* th[kMaxPar];
    int th_scalar[kMaxPar];
    lf7_inv12_t inv12 = nullptr;
    lf7_clamp_t clamp = nullptr;
    d7_sc_t sc = nullptr;

    inline void step(double e_t, int row1, double* s_t, double* c_t) const {
        double h, h1, h2;
        inv12(link, e_t, &h, &h1, &h2);
        double thl[kMaxPar];
        for (int j = 0; j < npar; ++j) {
            thl[j] = th_scalar[j] ? th[j][0] : th[j][row1 - 1];
        }
        thl[k] = clamp(h, lwr, upr);
        double out[2];
        sc(fam, k, y[row1 - 1], thl, out);
        // the same composition, in the same order, as distrib_kernel():
        // g * h1, and h * h1 * h1 + g * h2
        *s_t = out[0] * h1;
        *c_t = out[1] * h1 * h1 + out[0] * h2;
    }
};

FastCtx fast_ctx(SEXP fastS) {
    FastCtx c;
    if (Rf_isNull(fastS)) return c;
    List fast(fastS);
    std::string fam = as<std::string>(fast["family"]);
    std::string lnk = as<std::string>(fast["link"]);
    d7_id_t d7id = (d7_id_t) R_GetCCallable("distributions7", "d7_scalar_id");
    lf7_id_t lfid = (lf7_id_t) R_GetCCallable("linkfunctions7",
                                              "lf7_scalar_id");
    c.fam = d7id(fam.c_str());
    c.link = lfid(lnk.c_str());
    List th = fast["theta"];
    if (c.fam < 0 || c.link < 0 || th.size() > kMaxPar) return c;
    c.k = as<int>(fast["k"]) - 1;
    NumericVector bounds = fast["bounds"];
    c.lwr = bounds[0];
    c.upr = bounds[1];
    NumericVector y = fast["y"];
    c.y = y.begin();
    c.npar = th.size();
    for (int j = 0; j < c.npar; ++j) {
        NumericVector v = th[j];
        c.th[j] = v.begin();
        c.th_scalar[j] = (v.size() == 1);
    }
    c.inv12 = (lf7_inv12_t) R_GetCCallable("linkfunctions7", "lf7_inv12");
    c.clamp = (lf7_clamp_t) R_GetCCallable("linkfunctions7", "lf7_clamp");
    c.sc = (d7_sc_t) R_GetCCallable("distributions7", "d7_score_curv");
    c.ok = true;
    return c;
}

template <typename Body>
struct GroupWorker : public RcppParallel::Worker {
    const Body& body;
    explicit GroupWorker(const Body& b) : body(b) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t g = begin; g < end; ++g) body(g);
    }
};

// one group of the scalar route, everything a worker touches held as raw
// pointers extracted on the main thread
struct Grp {
    const int* rows;
    int m;
};

// and of the general route
struct GrpSub {
    const int* rows;
    const int* act;
    int m, na;
    const double* dom;               // m x na, column-major
    std::vector<const double*> da;   // p of them
    std::vector<const double*> dbl;  // q of them
    double f0;
    const double* df0;
};

} // namespace

// [[Rcpp::export]]
List gas_filter_cpp(NumericVector eta, List order, int p, int q,
                    double omega, NumericVector a, NumericVector b,
                    NumericMatrix db, double f0, NumericVector df0,
                    IntegerVector i_a, int np,
                    Function score, Function curvature,
                    Nullable<List> fast = R_NilValue, int threads = 1) {
    int n = eta.size();
    NumericVector eta_out(n);
    NumericMatrix jac(n, np);
    // the curvature at each step, returned so the adjoint's reverse pass
    // can read it instead of evaluating the callback a second time
    NumericVector curv(n);
    FastCtx fc = fast_ctx(fast);

    int ng = order.size();
    std::vector<Grp> grps(ng);
    std::vector<IntegerVector> keep(ng);   // keeps the R vectors alive
    for (int g = 0; g < ng; ++g) {
        keep[g] = as<IntegerVector>(order[g]);
        grps[g].rows = keep[g].begin();
        grps[g].m = keep[g].size();
    }

    const double* eta_p = eta.begin();
    double* out_p = eta_out.begin();
    double* jac_p = jac.begin();
    double* cv_p = curv.begin();
    const double* a_p = a.begin();
    const double* b_p = b.begin();
    const double* db_p = db.begin();
    const int dbn = db.nrow();
    const double* df0_p = df0.begin();
    const int* ia_p = i_a.begin();

    // one group's recursion; `use_fast` decides how the score and the
    // curvature are read, and nothing else differs between the routes
    auto run_group = [&](std::size_t gi, bool use_fast) {
        const Grp& G = grps[gi];
        int m = G.m;
        std::vector<double> f(m, 0.0), s(m, 0.0);
        std::vector<double> df(m * np, 0.0), ds(m * np, 0.0);
        std::vector<double> dft(np, 0.0);

        for (int t = 0; t < m; ++t) {
            double ft = omega;
            std::fill(dft.begin(), dft.end(), 0.0);
            dft[0] = 1.0;

            for (int i = 1; i <= p; ++i) {
                double s_lag = (t - i >= 0) ? s[t - i] : 0.0;
                ft += a_p[i - 1] * s_lag;
                if (t - i >= 0) {
                    for (int k = 0; k < np; ++k) {
                        dft[k] += a_p[i - 1] * ds[(t - i) * np + k];
                    }
                }
                dft[ia_p[i - 1] - 1] += s_lag;
            }

            for (int j = 1; j <= q; ++j) {
                double f_lag = (t - j >= 0) ? f[t - j] : f0;
                ft += b_p[j - 1] * f_lag;
                for (int k = 0; k < np; ++k) {
                    double df_lag = (t - j >= 0) ? df[(t - j) * np + k]
                                                 : df0_p[k];
                    dft[k] += b_p[j - 1] * df_lag + db_p[k * dbn + (j - 1)] * f_lag;
                }
            }

            f[t] = ft;
            for (int k = 0; k < np; ++k) df[t * np + k] = dft[k];

            int row = G.rows[t];
            double e_t = eta_p[row - 1] + ft;
            double s_t, c_t;
            if (use_fast) {
                fc.step(e_t, row, &s_t, &c_t);
            } else {
                s_t = as<double>(score(e_t, row));
                c_t = as<double>(curvature(e_t, row));
            }
            s[t] = s_t;
            for (int k = 0; k < np; ++k) ds[t * np + k] = c_t * dft[k];

            out_p[row - 1] = e_t;
            cv_p[row - 1] = c_t;
            for (int k = 0; k < np; ++k) {
                jac_p[k * n + (row - 1)] = dft[k];
            }
        }
    };

    if (fc.ok && threads > 1 && ng >= kMinGroupsPar) {
        auto body = [&](std::size_t gi) { run_group(gi, true); };
        GroupWorker<decltype(body)> w(body);
        RcppParallel::parallelFor(0, static_cast<std::size_t>(ng), w);
    } else {
        for (int g = 0; g < ng; ++g) run_group(g, fc.ok);
    }

    return List::create(_["eta"] = eta_out, _["jacobian"] = jac,
                        _["curv"] = curv);
}

// The general recursion of the submodel route: every coefficient of the
// filter read per observation, so a parameter carrying a subformula varies
// with t and the scalar kernel above does not apply. The derivative rows
// arrive already restricted to each group's active coordinates by the
// R-side preparation (one vectorized scatter per parameter per group), so
// the loop here allocates nothing per step.
//
// Each element of `groups` carries: `rows` (1-based observation indices in
// time order), `act` (1-based positions of the group's active coordinates
// among the term's np), `dom` (m x na, the derivative row of the level at
// each observation), `da` (p matrices m x na, the same for each loading),
// `dbl` (q matrices m x na, the derivative row of each autoregressive
// coefficient, already chained through Levinson-Durbin), and `f0`/`df0`
// (the starting level and its derivative, from the group's first
// observation). `om`, `A` and `B` are the per-observation values of the
// level, the loadings and the autoregressive coefficients.
//
// [[Rcpp::export]]
List gas_filter_sub_cpp(NumericVector eta, List groups, int p, int q,
                        NumericVector om, NumericMatrix A, NumericMatrix B,
                        int np, Function score, Function curvature,
                        Nullable<List> fast = R_NilValue, int threads = 1) {
    int n = eta.size();
    NumericVector eta_out(n);
    NumericMatrix jac(n, np);
    NumericVector curv(n);
    FastCtx fc = fast_ctx(fast);

    int ng = groups.size();
    std::vector<GrpSub> grps(ng);
    std::vector<List> keep(ng);            // keeps the R objects alive
    for (int g = 0; g < ng; ++g) {
        keep[g] = as<List>(groups[g]);
        List grp = keep[g];
        IntegerVector rows = grp["rows"];
        IntegerVector act = grp["act"];
        NumericMatrix dom = grp["dom"];
        List da = grp["da"];
        List dbl = grp["dbl"];
        GrpSub& G = grps[g];
        G.rows = rows.begin();
        G.act = act.begin();
        G.m = rows.size();
        G.na = act.size();
        G.dom = dom.begin();
        G.f0 = as<double>(grp["f0"]);
        NumericVector df0 = grp["df0"];
        G.df0 = df0.begin();
        for (int i = 0; i < p; ++i) {
            NumericMatrix dai = da[i];
            G.da.push_back(dai.begin());
        }
        for (int j = 0; j < q; ++j) {
            NumericMatrix dbj = dbl[j];
            G.dbl.push_back(dbj.begin());
        }
    }

    const double* eta_p = eta.begin();
    double* out_p = eta_out.begin();
    double* jac_p = jac.begin();
    double* cv_p = curv.begin();
    const double* om_p = om.begin();
    const double* A_p = A.begin();
    const double* B_p = B.begin();
    int An = A.nrow(), Bn = B.nrow();

    auto run_group = [&](std::size_t gi, bool use_fast) {
        const GrpSub& G = grps[gi];
        int m = G.m, na = G.na;
        std::vector<double> f(m, 0.0), s(m, 0.0);
        std::vector<double> df(m * na, 0.0), ds(m * na, 0.0);
        std::vector<double> dft(na, 0.0);

        for (int t = 0; t < m; ++t) {
            int r = G.rows[t];
            double ft = om_p[r - 1];
            for (int k = 0; k < na; ++k) dft[k] = G.dom[k * m + t];

            for (int i = 1; i <= p; ++i) {
                double a_r = A_p[(i - 1) * An + (r - 1)];
                const double* dai = G.da[i - 1];
                double s_lag = (t - i >= 0) ? s[t - i] : 0.0;
                ft += a_r * s_lag;
                if (t - i >= 0) {
                    for (int k = 0; k < na; ++k) {
                        dft[k] += a_r * ds[(t - i) * na + k];
                    }
                }
                for (int k = 0; k < na; ++k) {
                    dft[k] += s_lag * dai[k * m + t];
                }
            }

            for (int j = 1; j <= q; ++j) {
                double b_r = B_p[(j - 1) * Bn + (r - 1)];
                const double* dbj = G.dbl[j - 1];
                double f_lag = (t - j >= 0) ? f[t - j] : G.f0;
                ft += b_r * f_lag;
                for (int k = 0; k < na; ++k) {
                    double df_lag = (t - j >= 0) ? df[(t - j) * na + k]
                                                 : G.df0[k];
                    dft[k] += b_r * df_lag + f_lag * dbj[k * m + t];
                }
            }

            f[t] = ft;
            for (int k = 0; k < na; ++k) df[t * na + k] = dft[k];

            double e_t = eta_p[r - 1] + ft;
            double s_t, c_t;
            if (use_fast) {
                fc.step(e_t, r, &s_t, &c_t);
            } else {
                s_t = as<double>(score(e_t, r));
                c_t = as<double>(curvature(e_t, r));
            }
            s[t] = s_t;
            for (int k = 0; k < na; ++k) ds[t * na + k] = c_t * dft[k];

            out_p[r - 1] = e_t;
            cv_p[r - 1] = c_t;
            for (int k = 0; k < na; ++k) {
                jac_p[(G.act[k] - 1) * n + (r - 1)] = dft[k];
            }
        }
    };

    if (fc.ok && threads > 1 && ng >= kMinGroupsPar) {
        auto body = [&](std::size_t gi) { run_group(gi, true); };
        GroupWorker<decltype(body)> w(body);
        RcppParallel::parallelFor(0, static_cast<std::size_t>(ng), w);
    } else {
        for (int g = 0; g < ng; ++g) run_group(g, fc.ok);
    }

    return List::create(_["eta"] = eta_out, _["jacobian"] = jac,
                        _["curv"] = curv);
}
