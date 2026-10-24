#include <Rcpp.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "gas_par.h"
using namespace Rcpp;

// The joint score-driven recursion with coefficients that vary by
// observation, and its reverse pass, compiled. It is the per-row analogue of
// gas_joint_filter_cpp() and gas_joint_adjoint_cpp() (gas_joint.cpp): the
// family is read through distributions7's joint registry and the links
// through linkfunctions7's scalar entries, the driving vector and its
// derivative are computed as there, and the derivative of the levels is
// carried, as in gas_filter_sub_cpp(), on each group's active columns of the
// joint parameter vector instead of on all of them. The level, the loadings
// and the autoregressive coefficients are read at the row whose equation they
// enter; the cross loadings are scalars.
//
// The context (family, links, scaling) is the one of gas_joint_ctx.h,
// shared with gas_joint.cpp.
//
// Every index arriving from R is 1-based.

#include "gas_joint_ctx.h"
using namespace gasj;

namespace {

// The per-row coefficients of one member, as plain vectors so that a worker
// reads no R object: om (n), A (n x p, column-major), B (n x q), the partners'
// positions (1-based) and the scalar cross loadings (npk x p).
struct SubCoef {
    int p = 0, q = 0, n = 0;
    std::vector<double> om, A, B, cross;
    std::vector<int> partners;
};

std::vector<SubCoef> sub_coefs(const List& mem, bool need_om) {
    std::vector<SubCoef> out(mem.size());
    for (int k = 0; k < mem.size(); ++k) {
        List c = mem[k];
        SubCoef& J = out[k];
        J.p = as<int>(c["p"]);
        J.q = as<int>(c["q"]);
        if (need_om) {
            J.om = as<std::vector<double> >(c["om"]);
            J.n = J.om.size();
        }
        NumericMatrix A = as<NumericMatrix>(c["A"]);
        NumericMatrix B = as<NumericMatrix>(c["B"]);
        J.A.assign(A.begin(), A.end());
        J.B.assign(B.begin(), B.end());
        J.partners = as<std::vector<int> >(c["partners"]);
        NumericMatrix cr = as<NumericMatrix>(c["cross"]);
        J.cross.assign(cr.begin(), cr.end());
    }
    return out;
}

// one member of one group: the derivative rows of the level, of each loading
// and of each autoregressive coefficient on the group's active columns
// (m x na, column-major), the starting level and its derivative, and the
// positions among the active columns of the cross loadings (npk x p)
struct SubMem {
    const double* dom;
    std::vector<const double*> da, dbl;
    double f0;
    const double* df0;
    std::vector<int> ccol;
};

struct SubGrp {
    const int* rows;
    const int* act;
    int m, na;
    std::vector<SubMem> mem;
};

} // namespace

// The joint filter with per-row coefficients on the compiled route, or NULL
// where the context is inert (where gas_joint_filter_cpp() returns NULL).
// eta: n x K, the static predictors of the filtered parameters; mem: one
// list per member with p, q, om, A, B, partners, cross; groups: one list per
// group with rows, act (the group's active columns of the joint vector) and
// mem (one list per member with dom, da, dbl, f0, df0, ccol); M: the length
// of the joint parameter vector.
// [[Rcpp::export]]
SEXP gas_joint_sub_filter_cpp(NumericMatrix eta, List mem, List groups, int M,
                              Nullable<List> fast = R_NilValue,
                              int threads = 1) {
    JointCtx ctx = joint_ctx(fast.isNotNull() ? SEXP(fast.get()) : R_NilValue,
                         eta);
    if (!ctx.ok) return R_NilValue;
    const int n = eta.nrow(), K = eta.ncol();
    std::vector<SubCoef> co = sub_coefs(mem, true);
    NumericMatrix eta_out(n, K);
    std::vector<NumericMatrix> jac;
    for (int k = 0; k < K; ++k) jac.push_back(NumericMatrix(n, M));
    std::vector<double*> jac_p(K);
    for (int k = 0; k < K; ++k) jac_p[k] = jac[k].begin();
    NumericVector cv(n * K * K);
    double* cv_p = cv.begin();
    double* out_p = eta_out.begin();
    const double* eg = eta.begin();

    const int ng = groups.size();
    std::vector<SubGrp> grps(ng);
    // keeps the R objects alive while the workers hold their pointers
    std::vector<IntegerVector> keep_i;
    std::vector<NumericVector> keep_v;
    std::vector<NumericMatrix> keep_m;
    keep_i.reserve(2 * ng + 4 * ng * K);
    keep_v.reserve(2 * ng * K + 4);
    keep_m.reserve(static_cast<std::size_t>(ng) * K * 8 + 8);
    for (int g = 0; g < ng; ++g) {
        List grp = groups[g];
        SubGrp& G = grps[g];
        IntegerVector rows = grp["rows"];
        IntegerVector act = grp["act"];
        keep_i.push_back(rows);
        keep_i.push_back(act);
        G.rows = rows.begin();
        G.act = act.begin();
        G.m = rows.size();
        G.na = act.size();
        List gm = grp["mem"];
        G.mem.resize(K);
        for (int k = 0; k < K; ++k) {
            List mk = gm[k];
            SubMem& S = G.mem[k];
            NumericMatrix dom = mk["dom"];
            keep_m.push_back(dom);
            S.dom = dom.begin();
            List da = mk["da"];
            List dbl = mk["dbl"];
            for (int i = 0; i < co[k].p; ++i) {
                NumericMatrix x = da[i];
                keep_m.push_back(x);
                S.da.push_back(x.begin());
            }
            for (int j = 0; j < co[k].q; ++j) {
                NumericMatrix x = dbl[j];
                keep_m.push_back(x);
                S.dbl.push_back(x.begin());
            }
            S.f0 = as<double>(mk["f0"]);
            NumericVector df0 = mk["df0"];
            keep_v.push_back(df0);
            S.df0 = df0.begin();
            S.ccol = as<std::vector<int> >(mk["ccol"]);
        }
    }
    std::vector<int> bad(ng, 0);

    auto run_group = [&](std::size_t gi) {
        const SubGrp& Gp = grps[gi];
        const int T = Gp.m, na = Gp.na;
        const int* rows = Gp.rows;
        const std::size_t NA = static_cast<std::size_t>(na);
        std::vector<double> f(K * T, 0.0), u(K * T, 0.0);
        std::vector<double> Fd(static_cast<std::size_t>(K) * T * NA, 0.0),
            Sd(static_cast<std::size_t>(K) * T * NA, 0.0);
        std::vector<double> Ft(NA);
        double ut[kJMaxK], C[kJMaxK * kJMaxK], ft_all[kJMaxK];
        for (int t = 0; t < T; ++t) {
            const int row = rows[t];
            const int r = row - 1;
            for (int k = 0; k < K; ++k) {
                const SubCoef& J = co[k];
                const SubMem& S = Gp.mem[k];
                double ft = J.om[r];
                for (int c = 0; c < na; ++c) Ft[c] = S.dom[c * T + t];
                const int npk = J.partners.size();
                for (int i = 1; i <= J.p; ++i) {
                    const double ai = J.A[static_cast<std::size_t>(i - 1) * n + r];
                    const double* dai = S.da[i - 1];
                    if (t - i < 0) continue;   // the lagged score is zero
                    const std::size_t lag = static_cast<std::size_t>(k + K * (t - i));
                    const double sl = u[lag];
                    ft += ai * sl;
                    const double* sd = &Sd[lag * NA];
                    for (int c = 0; c < na; ++c) {
                        Ft[c] += ai * sd[c] + sl * dai[c * T + t];
                    }
                    for (int j = 0; j < npk; ++j) {
                        const int l = J.partners[j] - 1;
                        const double cji = J.cross[j + npk * (i - 1)];
                        const std::size_t ll = static_cast<std::size_t>(l + K * (t - i));
                        ft += cji * u[ll];
                        const double* sdl = &Sd[ll * NA];
                        for (int c = 0; c < na; ++c) Ft[c] += cji * sdl[c];
                        Ft[S.ccol[j + npk * (i - 1)] - 1] += u[ll];
                    }
                }
                for (int j = 1; j <= J.q; ++j) {
                    const double bj = J.B[static_cast<std::size_t>(j - 1) * n + r];
                    const double* dbj = S.dbl[j - 1];
                    double fl;
                    const double* fd;
                    if (t - j >= 0) {
                        const std::size_t lag = static_cast<std::size_t>(k + K * (t - j));
                        fl = f[lag];
                        fd = &Fd[lag * NA];
                    } else {
                        fl = S.f0;
                        fd = S.df0;
                    }
                    ft += bj * fl;
                    for (int c = 0; c < na; ++c) {
                        Ft[c] += bj * fd[c] + fl * dbj[c * T + t];
                    }
                }
                const std::size_t at = static_cast<std::size_t>(k + K * t);
                f[at] = ft;
                ft_all[k] = ft;
                std::copy(Ft.begin(), Ft.end(), &Fd[at * NA]);
            }
            if (!ctx.step(row, ft_all, ut, C)) {
                bad[gi] = 1;
                return;
            }
            for (int l = 0; l < K; ++l) {
                u[l + K * t] = ut[l];
                for (int k = 0; k < K; ++k) {
                    cv_p[r + n * (l + K * k)] = C[l + K * k];
                }
            }
            for (int l = 0; l < K; ++l) {
                double* sd = &Sd[static_cast<std::size_t>(l + K * t) * NA];
                for (int c = 0; c < na; ++c) {
                    double v = 0.0;
                    for (int k = 0; k < K; ++k) {
                        v += C[l + K * k] * Fd[static_cast<std::size_t>(k + K * t) * NA + c];
                    }
                    sd[c] = v;
                }
            }
            for (int k = 0; k < K; ++k) {
                out_p[r + n * k] = eg[r + n * k] + f[k + K * t];
                const double* fd = &Fd[static_cast<std::size_t>(k + K * t) * NA];
                for (int c = 0; c < na; ++c) {
                    jac_p[k][r + static_cast<std::size_t>(n) * (Gp.act[c] - 1)] = fd[c];
                }
            }
        }
    };

    mt7::GroupWorker<decltype(run_group)> w(run_group);
    if (ctx.thread_safe && threads > 1 && ng >= mt7::kMinGroupsPar) {
        RcppParallel::parallelFor(0, static_cast<std::size_t>(ng), w, 1,
                                  threads);
    } else {
        w(0, static_cast<std::size_t>(ng));
    }
    for (int g = 0; g < ng; ++g) {
        if (bad[g]) {
            stop("the expected information of the filtered parameters is not positive definite.");
        }
    }
    cv.attr("dim") = IntegerVector::create(n, K, K);
    List J(K);
    for (int k = 0; k < K; ++k) J[k] = jac[k];
    return List::create(_["eta"] = eta_out, _["jacobian"] = J,
                        _["curv"] = cv);
}

// The reverse pass of the joint recursion with per-row coefficients, from the
// curvature that the filter read: the derivative of a caller's objective in
// the static predictors (deta) and in the driving scores (dscore), given its
// direct derivative g (n x K) in the filtered predictors. The coefficients
// are read at the row whose equation they enter, and no derivative of them
// is needed.
// [[Rcpp::export]]
List gas_joint_sub_adjoint_cpp(NumericVector curv, NumericMatrix g, List mem,
                               List order) {
    const int n = g.nrow(), K = g.ncol();
    std::vector<SubCoef> co = sub_coefs(mem, false);
    NumericMatrix deta(n, K), dscore(n, K);
    for (int gi = 0; gi < order.size(); ++gi) {
        IntegerVector rows = order[gi];
        const int T = rows.size();
        std::vector<double> fb(K * T, 0.0), ub(K * T, 0.0);
        for (int t = T - 1; t >= 0; --t) {
            const int r = rows[t] - 1;
            for (int k = 0; k < K; ++k) {
                double eb = g(r, k);
                for (int l = 0; l < K; ++l) eb += curv[r + n * (l + K * k)] * ub[l + K * t];
                fb[k + K * t] += eb;
                deta(r, k) = eb;
                dscore(r, k) = ub[k + K * t];
            }
            for (int k = 0; k < K; ++k) {
                const SubCoef& J = co[k];
                const double fbk = fb[k + K * t];
                for (int j = 1; j <= J.q; ++j) {
                    if (t - j >= 0) {
                        fb[k + K * (t - j)] += J.B[static_cast<std::size_t>(j - 1) * n + r] * fbk;
                    }
                }
                const int npk = J.partners.size();
                for (int i = 1; i <= J.p; ++i) {
                    if (t - i < 0) continue;
                    ub[k + K * (t - i)] += J.A[static_cast<std::size_t>(i - 1) * n + r] * fbk;
                    for (int j = 0; j < npk; ++j) {
                        const int l = J.partners[j] - 1;
                        ub[l + K * (t - i)] += J.cross[j + npk * (i - 1)] * fbk;
                    }
                }
            }
        }
    }
    return List::create(_["deta"] = deta, _["dscore"] = dscore);
}
