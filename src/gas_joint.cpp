#include <Rcpp.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include <string>
#include <vector>
#include "gas_par.h"
using namespace Rcpp;

// The joint score-driven recursion of term_filter.GasJoint and its reverse
// pass, compiled. The filter reads the family through distributions7's
// joint registry (the gradient, the Hessian, the expected information and
// its derivative at one observation, one entry per order) and the links
// through linkfunctions7's scalar entries, all resolved once per call with
// R_GetCCallable, so the loop does not touch the R API. A family or a link
// that the registries do not cover makes the context inert, and the caller
// runs the R recursion with its callbacks instead. The arithmetic of the
// driving vector is the R route's (statmodels7's joint_callbacks()): the
// chain rule onto the link scale, the scaling I^-d of the filtered block with
// the symmetric root for d = 1/2, and C = S H + (dS) g. The reverse pass
// needs no family: it reads the curvature that the filter stored.
//
// Every index arriving from R is 1-based.

#include "gas_joint_ctx.h"
using namespace gasj;

namespace {


// the coefficients of one term of the joint recursion, from
// .gas_joint_coefs(), held as plain vectors so that a worker reads no R
// object
struct JCoef {
    int p = 0, q = 0, col_omega = 0;
    double omega = 0.0, f0 = 0.0;
    std::vector<double> a, b, db, df0, cross;   // db: q x nown, cross: np x p
    std::vector<int> own, col_a, partners, col_cross;
    int nown = 0;
};

std::vector<JCoef> joint_coefs(const List& cf) {
    std::vector<JCoef> out(cf.size());
    for (int k = 0; k < cf.size(); ++k) {
        List c = cf[k];
        JCoef& J = out[k];
        J.p = as<int>(c["p"]);
        J.q = as<int>(c["q"]);
        J.omega = as<double>(c["omega"]);
        J.f0 = as<double>(c["f0"]);
        J.col_omega = as<int>(c["col_omega"]);
        J.a = as<std::vector<double> >(c["a"]);
        J.b = as<std::vector<double> >(c["b"]);
        J.df0 = as<std::vector<double> >(c["df0"]);
        J.own = as<std::vector<int> >(c["own"]);
        J.col_a = as<std::vector<int> >(c["col_a"]);
        J.partners = as<std::vector<int> >(c["partners"]);
        J.nown = J.own.size();
        NumericMatrix db = as<NumericMatrix>(c["db"]);
        J.db.assign(db.begin(), db.end());
        NumericMatrix cr = as<NumericMatrix>(c["cross"]);
        J.cross.assign(cr.begin(), cr.end());
        IntegerMatrix cc = as<IntegerMatrix>(c["col_cross"]);
        J.col_cross.assign(cc.begin(), cc.end());
    }
    return out;
}

} // namespace

// The joint filter on the compiled route, or NULL where the context is inert.
// eta: n x K, the static predictors of the filtered parameters; cf:
// .gas_joint_coefs(); M: the length of the joint parameter vector; order:
// the groups' rows, 1-based, in time order.
// [[Rcpp::export]]
SEXP gas_joint_filter_cpp(NumericMatrix eta, List cf, int M, List order,
                          Nullable<List> fast = R_NilValue, int threads = 1) {
    JointCtx ctx = joint_ctx(fast.isNotNull() ? SEXP(fast.get()) : R_NilValue,
                             eta);
    if (!ctx.ok) return R_NilValue;
    const int n = eta.nrow(), K = eta.ncol();
    std::vector<JCoef> co = joint_coefs(cf);
    NumericMatrix eta_out(n, K);
    std::vector<NumericMatrix> jac;
    for (int k = 0; k < K; ++k) jac.push_back(NumericMatrix(n, M));
    std::vector<double*> jac_p(K);
    for (int k = 0; k < K; ++k) jac_p[k] = jac[k].begin();
    NumericVector cv(n * K * K);
    double* cv_p = cv.begin();
    double* out_p = eta_out.begin();
    const double* eg = eta.begin();

    int ng = order.size();
    std::vector<IntegerVector> keep(ng);
    std::vector<const int*> rows_p(ng);
    std::vector<int> len(ng);
    for (int g = 0; g < ng; ++g) {
        keep[g] = as<IntegerVector>(order[g]);
        rows_p[g] = keep[g].begin();
        len[g] = keep[g].size();
    }
    std::vector<int> bad(ng, 0);

    auto run_group = [&](std::size_t gi) {
        const int* rows = rows_p[gi];
        const int T = len[gi];
        std::vector<double> f(K * T, 0.0), u(K * T, 0.0);
        std::vector<double> Fd(static_cast<std::size_t>(K) * T * M, 0.0),
            Sd(static_cast<std::size_t>(K) * T * M, 0.0);
        std::vector<double> Ft(M), Fl(M);
        double ut[kJMaxK], C[kJMaxK * kJMaxK], ft_all[kJMaxK];
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < K; ++k) {
                const JCoef& J = co[k];
                double ft = J.omega;
                std::fill(Ft.begin(), Ft.end(), 0.0);
                Ft[J.col_omega - 1] = 1.0;
                for (int i = 1; i <= J.p; ++i) {
                    if (t - i < 0) continue;
                    const double ai = J.a[i - 1];
                    const std::size_t lag = static_cast<std::size_t>(k + K * (t - i));
                    ft += ai * u[lag];
                    const double* sd = &Sd[lag * M];
                    for (int m = 0; m < M; ++m) Ft[m] += ai * sd[m];
                    Ft[J.col_a[i - 1] - 1] += u[lag];
                    const int npk = J.partners.size();
                    for (int j = 0; j < npk; ++j) {
                        const int l = J.partners[j] - 1;
                        const double cji = J.cross[j + npk * (i - 1)];
                        const std::size_t ll = static_cast<std::size_t>(l + K * (t - i));
                        ft += cji * u[ll];
                        const double* sdl = &Sd[ll * M];
                        for (int m = 0; m < M; ++m) Ft[m] += cji * sdl[m];
                        Ft[J.col_cross[j + npk * (i - 1)] - 1] += u[ll];
                    }
                }
                for (int j = 1; j <= J.q; ++j) {
                    double fl;
                    if (t - j >= 0) {
                        const std::size_t lag = static_cast<std::size_t>(k + K * (t - j));
                        fl = f[lag];
                        const double* fd = &Fd[lag * M];
                        for (int m = 0; m < M; ++m) Fl[m] = fd[m];
                    } else {
                        fl = J.f0;
                        std::fill(Fl.begin(), Fl.end(), 0.0);
                        for (int o = 0; o < J.nown; ++o) Fl[J.own[o] - 1] = J.df0[o];
                    }
                    const double bj = J.b[j - 1];
                    ft += bj * fl;
                    for (int m = 0; m < M; ++m) Ft[m] += bj * Fl[m];
                    for (int o = 0; o < J.nown; ++o) {
                        Ft[J.own[o] - 1] += J.db[(j - 1) + J.q * o] * fl;
                    }
                }
                const std::size_t at = static_cast<std::size_t>(k + K * t);
                f[at] = ft;
                ft_all[k] = ft;
                std::copy(Ft.begin(), Ft.end(), &Fd[at * M]);
            }
            const int row = rows[t];
            if (!ctx.step(row, ft_all, ut, C)) {
                bad[gi] = 1;
                return;
            }
            for (int l = 0; l < K; ++l) {
                u[l + K * t] = ut[l];
                for (int k = 0; k < K; ++k) {
                    cv_p[(row - 1) + n * (l + K * k)] = C[l + K * k];
                }
            }
            for (int l = 0; l < K; ++l) {
                double* sd = &Sd[static_cast<std::size_t>(l + K * t) * M];
                for (int m = 0; m < M; ++m) {
                    double v = 0.0;
                    for (int k = 0; k < K; ++k) {
                        v += C[l + K * k] * Fd[static_cast<std::size_t>(k + K * t) * M + m];
                    }
                    sd[m] = v;
                }
            }
            for (int k = 0; k < K; ++k) {
                out_p[(row - 1) + n * k] = eg[(row - 1) + n * k] + f[k + K * t];
                const double* fd = &Fd[static_cast<std::size_t>(k + K * t) * M];
                for (int m = 0; m < M; ++m) jac_p[k][(row - 1) + static_cast<std::size_t>(n) * m] = fd[m];
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

// The reverse pass of the joint recursion from the curvature that the
// filter read: the derivative of a caller's objective in the static
// predictors (deta) and in the driving scores (dscore), given its direct
// derivative g (n x K) in the filtered predictors.
// [[Rcpp::export]]
List gas_joint_adjoint_cpp(NumericVector curv, NumericMatrix g, List cf,
                           List order) {
    const int n = g.nrow(), K = g.ncol();
    std::vector<JCoef> co = joint_coefs(cf);
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
                const JCoef& J = co[k];
                const double fbk = fb[k + K * t];
                for (int j = 1; j <= J.q; ++j) {
                    if (t - j >= 0) fb[k + K * (t - j)] += J.b[j - 1] * fbk;
                }
                const int npk = J.partners.size();
                for (int i = 1; i <= J.p; ++i) {
                    if (t - i < 0) continue;
                    ub[k + K * (t - i)] += J.a[i - 1] * fbk;
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
