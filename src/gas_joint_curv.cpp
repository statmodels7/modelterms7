#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "gas_par.h"
using namespace Rcpp;

// The second order of the joint score-driven filter, compiled: what
// .gas_joint_core() in R/gas_joint_core.R returns when it has no direction,
// the Jacobian of every filtered predictor in the unknowns and the second
// derivative of the filtered predictors contracted against the weights g.
//
// Each level and each driving score is carried as a value, a gradient (m) and
// a Hessian (m x m) in the unknowns, as the R code carries them as jets, and
// the arrays are written out per order: the product of two quantities is
// Leibniz's rule at orders zero, one and two, and the driving score's
// derivatives come from the chain rule for a function of the predictors of
// every equation,
//   d1 = Vm U1,   d2 = Vm U2 Vm' + sum_a U1[a] E_a'' ,
// with Vm the gradients of the predictors side by side and E_a'' the Hessian
// of predictor a (non-zero only for a filtered equation). The charts of the
// terms (the coefficients and their first and second derivatives in the
// unknowns) are prepared in R by .gas_chart_prep() and read here as plain
// vectors, so a worker touches no R object.
//
// Each group of rows keeps the last L + 1 time points of every level and
// driving score in a ring (L the longest lag of any term), so the memory of
// a group does not grow with its length. The groups run over threads in
// blocks; the curvature of each group is stored apart and summed in group
// order, so the result does not depend on the thread count.
//
// Every index arriving from R is 1-based.

namespace {

// the chart of one term lifted to the m unknowns, from .gas_chart_prep():
// omega and the starting level f0 with their gradients (m) and Hessians
// (m x m, column-major); a and b the coefficients with gradients stored as
// m x p and m x q, Hessians as m x m x p and m x m x q
struct JCurvChart {
    double omega = 0.0, f0 = 0.0;
    std::vector<double> om_u, om_uu, f0_u, f0_uu;
    std::vector<double> a, a_u, a_uu, b, b_u, b_uu;
};

// the coefficients of one term of the joint recursion, from
// .gas_joint_coefs(): the orders, the partners (cross loadings) and the
// columns at which a cross loading sits
struct JCurvCoef {
    int p = 0, q = 0;
    std::vector<int> partners, col_cross;   // col_cross: np x p
    std::vector<double> cross;              // np x p
};

// J += A B for two quantities with value, gradient and Hessian: Leibniz's
// rule through the second order. A2 may be null (a quantity that is linear in
// the unknowns has no Hessian).
inline void prod_add(int m, double Av, const double* A1, const double* A2,
                     double Bv, const double* B1, const double* B2,
                     double* Jv, double* J1, double* J2) {
    *Jv += Av * Bv;
    for (int i = 0; i < m; ++i) J1[i] += Av * B1[i] + Bv * A1[i];
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            double v = Av * B2[i + m * j] + A1[i] * B1[j] + B1[i] * A1[j];
            if (A2) v += Bv * A2[i + m * j];
            J2[i + m * j] += v;
        }
    }
}

} // namespace

// [[Rcpp::export]]
List gas_joint_curv_cpp(List seed, List charts, List cf, IntegerVector fp,
                        NumericMatrix g, NumericMatrix u_at,
                        NumericVector U1, NumericVector U2, List order,
                        int m, int zoff, int threads = 1) {
    const int K = g.ncol(), n = g.nrow();
    const int P = seed.size();
    if (charts.size() != K || cf.size() != K || fp.size() != K ||
        u_at.nrow() != n || u_at.ncol() != K) {
        stop("the inputs of the compiled joint curvature do not agree in size.");
    }
    if (U1.size() != static_cast<R_xlen_t>(n) * K * P ||
        U2.size() != static_cast<R_xlen_t>(n) * K * P * P) {
        stop("'U1' and 'U2' must have dimensions c(n, K, P) and c(n, K, P, P).");
    }
    const std::size_t mm = static_cast<std::size_t>(m) * m;
    std::vector<NumericMatrix> seedm(P);
    std::vector<const double*> seed_p(P);
    for (int a = 0; a < P; ++a) {
        seedm[a] = as<NumericMatrix>(seed[a]);
        if (seedm[a].nrow() != n || seedm[a].ncol() < m) {
            stop("each matrix of 'seed' must have one row per observation and one column per unknown.");
        }
        seed_p[a] = seedm[a].begin();
    }
    std::vector<int> fpe(K);
    for (int k = 0; k < K; ++k) {
        fpe[k] = fp[k] - 1;
        if (fpe[k] < 0 || fpe[k] >= P) {
            stop("'fp' must give the equation of each filtered parameter.");
        }
    }

    std::vector<JCurvChart> ch(K);
    std::vector<JCurvCoef> co(K);
    int L = 0;
    for (int k = 0; k < K; ++k) {
        List c = charts[k];
        List f = cf[k];
        JCurvChart& C = ch[k];
        JCurvCoef& J = co[k];
        J.p = as<int>(f["p"]);
        J.q = as<int>(f["q"]);
        J.partners = as<std::vector<int> >(f["partners"]);
        const std::size_t np = J.partners.size();
        J.cross = as<std::vector<double> >(f["cross"]);
        J.col_cross = as<std::vector<int> >(f["col_cross"]);
        C.omega = as<double>(c["omega"]);
        C.f0 = as<double>(c["f0"]);
        C.om_u = as<std::vector<double> >(c["om_u"]);
        C.om_uu = as<std::vector<double> >(c["om_uu"]);
        C.f0_u = as<std::vector<double> >(c["f0_u"]);
        C.f0_uu = as<std::vector<double> >(c["f0_uu"]);
        C.a = as<std::vector<double> >(c["a"]);
        C.a_u = as<std::vector<double> >(c["a_u"]);
        C.a_uu = as<std::vector<double> >(c["a_uu"]);
        C.b = as<std::vector<double> >(c["b"]);
        C.b_u = as<std::vector<double> >(c["b_u"]);
        C.b_uu = as<std::vector<double> >(c["b_uu"]);
        const std::size_t p = J.p, q = J.q;
        if (C.om_u.size() != static_cast<std::size_t>(m) || C.om_uu.size() != mm ||
            C.f0_u.size() != static_cast<std::size_t>(m) || C.f0_uu.size() != mm ||
            C.a.size() != p || C.a_u.size() != p * m || C.a_uu.size() != p * mm ||
            C.b.size() != q || C.b_u.size() != q * m || C.b_uu.size() != q * mm ||
            J.cross.size() != np * p || J.col_cross.size() != np * p) {
            stop("the chart of a term does not agree with its orders and the number of unknowns.");
        }
        for (std::size_t j = 0; j < np; ++j) {
            if (J.partners[j] < 1 || J.partners[j] > K) {
                stop("a partner of a filtered parameter is out of range.");
            }
        }
        for (std::size_t i = 0; i < J.col_cross.size(); ++i) {
            if (J.col_cross[i] < 1 || zoff + J.col_cross[i] > m) {
                stop("a cross loading sits outside the unknowns.");
            }
        }
        L = std::max(L, std::max(J.p, J.q));
    }
    const int D = L + 1;   // the ring's depth

    NumericMatrix curvature(m, m);
    std::vector<NumericMatrix> jac;
    for (int k = 0; k < K; ++k) jac.push_back(NumericMatrix(n, m));
    std::vector<double*> jac_p(K);
    for (int k = 0; k < K; ++k) jac_p[k] = jac[k].begin();
    const double* gp = g.begin();
    const double* up = u_at.begin();
    const double* U1p = U1.begin();
    const double* U2p = U2.begin();

    const int ng = order.size();
    std::vector<IntegerVector> keep(ng);
    std::vector<const int*> rows_p(ng);
    std::vector<int> len(ng);
    for (int gi = 0; gi < ng; ++gi) {
        keep[gi] = as<IntegerVector>(order[gi]);
        rows_p[gi] = keep[gi].begin();
        len[gi] = keep[gi].size();
        for (int t = 0; t < len[gi]; ++t) {
            if (rows_p[gi][t] < 1 || rows_p[gi][t] > n) {
                stop("'order' holds a row outside the data.");
            }
        }
    }

    const int blockg = 256;
    std::vector<double> Wblk;
    std::size_t b0 = 0;

    auto run_group = [&](std::size_t gi) {
        const int* rows = rows_p[gi];
        const int T = len[gi];
        double* Wg = &Wblk[(gi - b0) * mm];
        // the rings: level (f) and driving score (u), one slot per
        // (parameter, time mod D)
        const std::size_t ns = static_cast<std::size_t>(K) * D;
        std::vector<double> fv(ns, 0.0), uv(ns, 0.0);
        std::vector<double> fd1(ns * m, 0.0), ud1(ns * m, 0.0);
        std::vector<double> fd2(ns * mm, 0.0), ud2(ns * mm, 0.0);
        std::vector<double> ecol(m, 0.0);
        std::vector<double> Ed1(static_cast<std::size_t>(P) * m);
        std::vector<const double*> Ed2(P);
        std::vector<double> TM(static_cast<std::size_t>(m) * P);
        for (int t = 0; t < T; ++t) {
            const int ts = t % D;
            // the levels of every filtered parameter at this time
            for (int k = 0; k < K; ++k) {
                const JCurvChart& C = ch[k];
                const JCurvCoef& J = co[k];
                const std::size_t cur = static_cast<std::size_t>(k + K * ts);
                double Jv = C.omega;
                double* J1 = &fd1[cur * m];
                double* J2 = &fd2[cur * mm];
                std::copy(C.om_u.begin(), C.om_u.end(), J1);
                std::copy(C.om_uu.begin(), C.om_uu.end(), J2);
                const int np = static_cast<int>(J.partners.size());
                for (int i = 1; i <= J.p; ++i) {
                    if (t - i < 0) continue;
                    const int ls = (t - i) % D;
                    // the score loading times the parameter's own lagged score
                    {
                        const std::size_t lag = static_cast<std::size_t>(k + K * ls);
                        prod_add(m, C.a[i - 1], &C.a_u[static_cast<std::size_t>(i - 1) * m],
                                 &C.a_uu[static_cast<std::size_t>(i - 1) * mm],
                                 uv[lag], &ud1[lag * m], &ud2[lag * mm],
                                 &Jv, J1, J2);
                    }
                    // the cross loadings times the partners' lagged scores; a
                    // cross loading rides the identity link, so its gradient
                    // is a one at its own column and it has no Hessian
                    for (int j = 0; j < np; ++j) {
                        const int l = J.partners[j] - 1;
                        const std::size_t lag = static_cast<std::size_t>(l + K * ls);
                        const std::size_t ix = static_cast<std::size_t>(j) + static_cast<std::size_t>(np) * (i - 1);
                        const int col = zoff + J.col_cross[ix] - 1;
                        ecol[col] = 1.0;
                        prod_add(m, J.cross[ix], ecol.data(), nullptr,
                                 uv[lag], &ud1[lag * m], &ud2[lag * mm],
                                 &Jv, J1, J2);
                        ecol[col] = 0.0;
                    }
                }
                for (int j = 1; j <= J.q; ++j) {
                    double lv;
                    const double* l1;
                    const double* l2;
                    if (t - j >= 0) {
                        const std::size_t lag = static_cast<std::size_t>(k + K * ((t - j) % D));
                        lv = fv[lag];
                        l1 = &fd1[lag * m];
                        l2 = &fd2[lag * mm];
                    } else {
                        lv = C.f0;
                        l1 = C.f0_u.data();
                        l2 = C.f0_uu.data();
                    }
                    prod_add(m, C.b[j - 1], &C.b_u[static_cast<std::size_t>(j - 1) * m],
                             &C.b_uu[static_cast<std::size_t>(j - 1) * mm],
                             lv, l1, l2, &Jv, J1, J2);
                }
                fv[cur] = Jv;
            }
            const int r = rows[t] - 1;
            // the predictors of every equation at this row: the static rows,
            // and for a filtered equation the level's gradient on top and
            // its Hessian
            for (int a = 0; a < P; ++a) {
                for (int i = 0; i < m; ++i) {
                    Ed1[static_cast<std::size_t>(a) * m + i] =
                        seed_p[a][r + static_cast<std::size_t>(n) * i];
                }
                Ed2[a] = nullptr;
            }
            for (int k = 0; k < K; ++k) {
                const int a = fpe[k];
                const std::size_t cur = static_cast<std::size_t>(k + K * ts);
                const double* J1 = &fd1[cur * m];
                const double* J2 = &fd2[cur * mm];
                double* e1 = &Ed1[static_cast<std::size_t>(a) * m];
                for (int i = 0; i < m; ++i) e1[i] += J1[i];
                Ed2[a] = J2;
                for (int i = 0; i < m; ++i) {
                    jac_p[k][r + static_cast<std::size_t>(n) * i] = e1[i];
                }
                const double gk = gp[r + static_cast<std::size_t>(n) * k];
                for (std::size_t i = 0; i < mm; ++i) Wg[i] += gk * J2[i];
            }
            // the driving score of every filtered parameter at this row
            for (int l = 0; l < K; ++l) {
                const std::size_t cur = static_cast<std::size_t>(l + K * ts);
                const double* u1 = U1p + r + static_cast<std::size_t>(n) * l;
                const double* u2 = U2p + r + static_cast<std::size_t>(n) * l;
                const std::size_t s1 = static_cast<std::size_t>(n) * K;
                uv[cur] = up[r + static_cast<std::size_t>(n) * l];
                double* d1 = &ud1[cur * m];
                double* d2 = &ud2[cur * mm];
                // d1 = Vm U1
                for (int i = 0; i < m; ++i) {
                    double v = 0.0;
                    for (int a = 0; a < P; ++a) {
                        v += Ed1[static_cast<std::size_t>(a) * m + i] * u1[s1 * a];
                    }
                    d1[i] = v;
                }
                // d2 = (Vm U2) Vm' + sum_a U1[a] E_a''
                for (int b = 0; b < P; ++b) {
                    for (int i = 0; i < m; ++i) {
                        double v = 0.0;
                        for (int a = 0; a < P; ++a) {
                            v += Ed1[static_cast<std::size_t>(a) * m + i] *
                                 u2[s1 * (a + static_cast<std::size_t>(P) * b)];
                        }
                        TM[i + static_cast<std::size_t>(m) * b] = v;
                    }
                }
                for (int j = 0; j < m; ++j) {
                    for (int i = 0; i < m; ++i) {
                        double v = 0.0;
                        for (int b = 0; b < P; ++b) {
                            v += TM[i + static_cast<std::size_t>(m) * b] *
                                 Ed1[static_cast<std::size_t>(b) * m + j];
                        }
                        d2[i + m * j] = v;
                    }
                }
                for (int a = 0; a < P; ++a) {
                    if (!Ed2[a]) continue;
                    const double ua = u1[s1 * a];
                    for (std::size_t i = 0; i < mm; ++i) d2[i] += ua * Ed2[a][i];
                }
            }
        }
    };

    mt7::GroupWorker<decltype(run_group)> w(run_group);
    std::vector<double> W(mm, 0.0);
    for (b0 = 0; b0 < static_cast<std::size_t>(ng); b0 += blockg) {
        const std::size_t b1 = std::min(static_cast<std::size_t>(ng),
                                        b0 + static_cast<std::size_t>(blockg));
        Wblk.assign((b1 - b0) * mm, 0.0);
        if (threads > 1 && ng >= mt7::kMinGroupsPar) {
            RcppParallel::parallelFor(b0, b1, w, 1, threads);
        } else {
            w(b0, b1);
        }
        for (std::size_t gi = 0; gi < b1 - b0; ++gi) {
            for (std::size_t i = 0; i < mm; ++i) W[i] += Wblk[gi * mm + i];
        }
    }
    // the entries and their transposes collect the same terms in a different
    // order, so the sum is made symmetric exactly
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            curvature(i, j) = 0.5 * (W[i + m * j] + W[j + m * i]);
        }
    }
    List J(K);
    for (int k = 0; k < K; ++k) J[k] = jac[k];
    return List::create(_["jacobian"] = J, _["curvature"] = curvature);
}
