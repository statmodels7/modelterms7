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
    // the third derivatives along each direction (nd blocks of m x m, and
    // for a and b nd x p and nd x q blocks) and, with two directions, the
    // fourth derivatives contracted against both; empty at the second order
    std::vector<double> om_t, a_t, b_t, f0_3, om_q, a_q, b_q, f0_4;
};

// the coefficients of one term of the joint recursion, from
// .gas_joint_coefs(): the orders, the partners (cross loadings) and the
// columns at which a cross loading sits
struct JCurvCoef {
    int p = 0, q = 0;
    // one entry per cross pair (partner xpar, lag xlag), in the order of
    // .gas_cross_params(); a pair's lag may exceed p
    std::vector<int> partners, col_cross, xpar, xlag;
    std::vector<double> cross;
    int nx = 0, xmax = 0;                   // the pairs, and the largest lag
    // a static start: the first level of each group is d0[group], with the
    // derivative minus the static row of its equation at the first row
    bool direct = false;
    std::vector<double> d0;
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

// reads and checks the cross pairs of one term from .gas_joint_coefs():
// xpar, xlag, cross and col_cross, one entry per pair
inline void read_cross_pairs(JCurvCoef& J, const List& f, int K, int zoff,
                             int m) {
    J.partners = as<std::vector<int> >(f["partners"]);
    J.xpar = as<std::vector<int> >(f["xpar"]);
    J.xlag = as<std::vector<int> >(f["xlag"]);
    J.cross = as<std::vector<double> >(f["cross"]);
    J.col_cross = as<std::vector<int> >(f["col_cross"]);
    J.nx = static_cast<int>(J.xpar.size());
    if (J.xlag.size() != J.xpar.size() || J.cross.size() != J.xpar.size() ||
        J.col_cross.size() != J.xpar.size()) {
        stop("the cross loadings of a term must have one partner, lag, value and column per pair.");
    }
    if (f.containsElementNamed("direct") && as<bool>(f["direct"])) {
        J.direct = true;
        J.d0 = as<std::vector<double> >(f["d0"]);
    }
    J.xmax = 0;
    for (int h = 0; h < J.nx; ++h) {
        if (J.xpar[h] < 1 || J.xpar[h] > K) {
            stop("a partner of a filtered parameter is out of range.");
        }
        if (J.xlag[h] < 1) {
            stop("the lag of a cross loading must be at least one.");
        }
        if (J.col_cross[h] < 1 || zoff + J.col_cross[h] > m) {
            stop("a cross loading sits outside the unknowns.");
        }
        J.xmax = std::max(J.xmax, J.xlag[h]);
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
        read_cross_pairs(J, f, K, zoff, m);
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
            C.b.size() != q || C.b_u.size() != q * m || C.b_uu.size() != q * mm) {
            stop("the chart of a term does not agree with its orders and the number of unknowns.");
        }
        // the rings must reach the largest cross lag as well
        L = std::max(L, std::max(std::max(J.p, J.q), J.xmax));
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
        // the static starts of the members: minus the static row of each
        // member's equation at the group's first row, and no Hessian
        std::vector<double> S1(static_cast<std::size_t>(K) * m, 0.0), Z2(mm, 0.0);
        for (int k = 0; k < K; ++k) {
            if (!co[k].direct || T == 0) continue;
            const int r1 = rows[0] - 1;
            for (int i = 0; i < m; ++i) {
                S1[static_cast<std::size_t>(k) * m + i] =
                    -seed_p[fpe[k]][r1 + static_cast<std::size_t>(n) * i];
            }
        }
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
                const double* s1 = &S1[static_cast<std::size_t>(k) * m];
                if (J.direct && t == 0) {
                    fv[cur] = J.d0[gi];
                    std::copy(s1, s1 + m, J1);
                    std::fill(J2, J2 + mm, 0.0);
                    continue;
                }
                std::copy(C.om_u.begin(), C.om_u.end(), J1);
                std::copy(C.om_uu.begin(), C.om_uu.end(), J2);
                for (int i = 1; i <= J.p; ++i) {
                    if (t - i < 0) continue;
                    const int ls = (t - i) % D;
                    // the score loading times the parameter's own lagged score
                    const std::size_t lag = static_cast<std::size_t>(k + K * ls);
                    prod_add(m, C.a[i - 1], &C.a_u[static_cast<std::size_t>(i - 1) * m],
                             &C.a_uu[static_cast<std::size_t>(i - 1) * mm],
                             uv[lag], &ud1[lag * m], &ud2[lag * mm],
                             &Jv, J1, J2);
                }
                // the cross loadings times the partners' lagged scores, one
                // term per pair at its own lag; a cross loading rides the
                // identity link, so its gradient is a one at its own column
                // and it has no Hessian
                for (int h = 0; h < J.nx; ++h) {
                    const int xl = J.xlag[h];
                    if (t - xl < 0) continue;
                    const std::size_t lag = static_cast<std::size_t>((J.xpar[h] - 1) + K * ((t - xl) % D));
                    const int col = zoff + J.col_cross[h] - 1;
                    ecol[col] = 1.0;
                    prod_add(m, J.cross[h], ecol.data(), nullptr,
                             uv[lag], &ud1[lag * m], &ud2[lag * mm],
                             &Jv, J1, J2);
                    ecol[col] = 0.0;
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
                    } else if (J.direct) {
                        lv = J.d0[gi];
                        l1 = s1;
                        l2 = Z2.data();
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

// ---------------------------------------------------------------------------
// The third and fourth orders of the joint filter, compiled: what
// .gas_joint_core() returns for one direction (term_third) and for two
// (term_fourth). Every level and driving score carries, besides the value, the
// gradient and the Hessian of the second order, its third derivative
// contracted against each direction (m x m) and, with two directions, its
// fourth derivative contracted against both (m x m), as the jets of the R
// recursion do. The arrays are written out per order: the product of two
// quantities is Leibniz's rule at order three (prod_add3) and at order four
// (prod_add4, the sixteen terms of .gas_prod4()), and the driving score's
// derivatives are the chain rule of .gas_jdrive() at each order. The
// contractions of a quantity that the products read (its gradient, Hessian and
// third derivative against the directions) are formed once, when the quantity
// is complete, and kept beside it in the ring.

namespace {

// a quantity seen through pointers: value, gradient, Hessian (null for zero),
// the third derivative along each direction (null for zero), the fourth
// (null for zero) and the contractions of jet_aux()
struct JetV {
    double v = 0.0;
    const double* d1 = nullptr;
    const double* d2 = nullptr;
    const double* d3[2] = {nullptr, nullptr};
    const double* d4 = nullptr;
    double ad[2] = {0.0, 0.0};                  // d1 . x_d
    const double* d2x[2] = {nullptr, nullptr};  // d2 x_d (m)
    double d2vw = 0.0;                          // x_0' d2 x_1
    const double* d3vw = nullptr;               // d3 along x_0, applied to x_1 (m)
};

// the contractions of a quantity against the directions x0 and x1 (the second
// one used only with two directions): d1 . x_d, d2 x_d, x0' d2 x1 and d3 x1
inline void jet_aux(int m, int nd, const double* x0, const double* x1,
                    const double* d1, const double* d2, const double* d3v,
                    double* ad, double* d2x, double* d2vw, double* d3vw) {
    const double* xs[2] = {x0, x1};
    for (int d = 0; d < nd; ++d) {
        const double* x = xs[d];
        double s = 0.0;
        for (int i = 0; i < m; ++i) s += d1[i] * x[i];
        ad[d] = s;
        for (int i = 0; i < m; ++i) {
            double t = 0.0;
            for (int j = 0; j < m; ++j) {
                t += d2[i + static_cast<std::size_t>(m) * j] * x[j];
            }
            d2x[static_cast<std::size_t>(d) * m + i] = t;
        }
    }
    if (nd == 2) {
        double s = 0.0;
        for (int i = 0; i < m; ++i) s += x0[i] * d2x[m + i];
        *d2vw = s;
        for (int i = 0; i < m; ++i) {
            double t = 0.0;
            for (int j = 0; j < m; ++j) {
                t += d3v[i + static_cast<std::size_t>(m) * j] * x1[j];
            }
            d3vw[i] = t;
        }
    }
}

// J3 += the third derivative of A B along each direction, J3 holding nd
// blocks of m x m. Leibniz's rule, as .gas_jprod() writes it. B is a lagged
// level or driving score or the starting level and has every derivative; A is
// a chart coefficient or a cross loading, whose missing derivatives are null.
inline void prod_add3(int m, int nd, const JetV& A, const JetV& B, double* J3) {
    const std::size_t mm = static_cast<std::size_t>(m) * m;
    for (int d = 0; d < nd; ++d) {
        double* J = J3 + static_cast<std::size_t>(d) * mm;
        const double* A1 = A.d1;
        const double* B1 = B.d1;
        const double* A2 = A.d2;
        const double* B2 = B.d2;
        const double* A3 = A.d3[d];
        const double* B3 = B.d3[d];
        const double* A2x = A.d2x[d];
        const double* B2x = B.d2x[d];
        const double ax = A.ad[d];
        const double bx = B.ad[d];
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < m; ++i) {
                const std::size_t ij = i + static_cast<std::size_t>(m) * j;
                double v = A.v * B3[ij] + ax * B2[ij] +
                           B2x[i] * A1[j] + A1[i] * B2x[j];
                if (A3) v += B.v * A3[ij];
                if (A2) v += bx * A2[ij];
                if (A2x) v += A2x[i] * B1[j] + B1[i] * A2x[j];
                J[ij] += v;
            }
        }
    }
}

// J4 += the fourth derivative of A B contracted against both directions: the
// sixteen terms of .gas_prod4(). The same conventions as prod_add3().
inline void prod_add4(int m, const JetV& A, const JetV& B, double* J4) {
    const double* A1 = A.d1;
    const double* B1 = B.d1;
    const double* A2 = A.d2;
    const double* B2 = B.d2;
    const double* A3v = A.d3[0];
    const double* A3w = A.d3[1];
    const double* B3v = B.d3[0];
    const double* B3w = B.d3[1];
    const double* A2v = A.d2x[0];
    const double* A2w = A.d2x[1];
    const double* B2v = B.d2x[0];
    const double* B2w = B.d2x[1];
    const double* A3vw = A.d3vw;
    const double* B3vw = B.d3vw;
    const double Aa = A.ad[0], Ab = A.ad[1], Ba = B.ad[0], Bb = B.ad[1];
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            const std::size_t ij = i + static_cast<std::size_t>(m) * j;
            double v = A.v * B.d4[ij] + Aa * B3w[ij] + Ab * B3v[ij] +
                       A.d2vw * B2[ij] +
                       A1[i] * B3vw[j] + B3vw[i] * A1[j];
            if (A2) v += B.d2vw * A2[ij];
            if (A2v) {
                v += A2v[i] * B2w[j] + B2w[i] * A2v[j] +
                     A2w[i] * B2v[j] + B2v[i] * A2w[j];
            }
            if (A3vw) v += A3vw[i] * B1[j] + B1[i] * A3vw[j];
            if (A3v) v += A3v[ij] * Bb;
            if (A3w) v += A3w[ij] * Ba;
            if (A.d4) v += A.d4[ij] * B.v;
            J4[ij] += v;
        }
    }
}

// out (m x P) = X (m x P) M (P x P)
inline void mat_mulp(int m, int P, double* out, const double* X, const double* M) {
    for (int b = 0; b < P; ++b) {
        for (int i = 0; i < m; ++i) {
            double s = 0.0;
            for (int a = 0; a < P; ++a) {
                s += X[i + static_cast<std::size_t>(m) * a] *
                     M[a + static_cast<std::size_t>(P) * b];
            }
            out[i + static_cast<std::size_t>(m) * b] = s;
        }
    }
}

// out (m x m) += X Y' for X and Y of size m x P
inline void add_xyt(int m, int P, double* out, const double* X, const double* Y) {
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            double s = 0.0;
            for (int b = 0; b < P; ++b) {
                s += X[i + static_cast<std::size_t>(m) * b] *
                     Y[j + static_cast<std::size_t>(m) * b];
            }
            out[i + static_cast<std::size_t>(m) * j] += s;
        }
    }
}

// the chart's coefficients seen as quantities, with their contractions
struct JChartJets {
    JetV f0;
    std::vector<JetV> a, b;
    std::vector<double> store;
};

// the ring of one kind of quantity (levels or driving scores): one slot per
// (parameter, time mod depth), with every derivative and the contractions
struct JRing {
    int m = 0, nd = 0;
    std::size_t mm = 0;
    std::vector<double> v, d1, d2, d3, d4, ad, d2x, d2vw, d3vw;
    void init(std::size_t ns, int m_, int nd_) {
        m = m_;
        nd = nd_;
        mm = static_cast<std::size_t>(m) * m;
        v.assign(ns, 0.0);
        d1.assign(ns * m, 0.0);
        d2.assign(ns * mm, 0.0);
        d3.assign(ns * nd * mm, 0.0);
        d4.assign(nd == 2 ? ns * mm : 0, 0.0);
        ad.assign(ns * nd, 0.0);
        d2x.assign(ns * nd * m, 0.0);
        d2vw.assign(nd == 2 ? ns : 0, 0.0);
        d3vw.assign(nd == 2 ? ns * m : 0, 0.0);
    }
    void aux(std::size_t s, const double* x0, const double* x1) {
        jet_aux(m, nd, x0, x1, &d1[s * m], &d2[s * mm], &d3[s * nd * mm],
                &ad[s * nd], &d2x[s * nd * m],
                nd == 2 ? &d2vw[s] : nullptr, nd == 2 ? &d3vw[s * m] : nullptr);
    }
    JetV view(std::size_t s) const {
        JetV V;
        V.v = v[s];
        V.d1 = &d1[s * m];
        V.d2 = &d2[s * mm];
        for (int d = 0; d < nd; ++d) {
            V.d3[d] = &d3[(s * nd + d) * mm];
            V.ad[d] = ad[s * nd + d];
            V.d2x[d] = &d2x[(s * nd + d) * m];
        }
        if (nd == 2) {
            V.d4 = &d4[s * mm];
            V.d2vw = d2vw[s];
            V.d3vw = &d3vw[s * m];
        }
        return V;
    }
};

// The recursion for nd = 1 or 2 directions. U3 and U4 are the driving
// scores' third and fourth derivatives in the predictors, arrays of dimension
// c(n, K, P, P, P) and c(n, K, P, P, P, P); U4 is read only with two
// directions.
List joint_run(int nd, List seed, List charts, List cf, IntegerVector fp,
               NumericMatrix g, NumericMatrix u_at, NumericVector U1,
               NumericVector U2, NumericVector U3, NumericVector U4,
               const std::vector<std::vector<double> >& dirs, List order,
               int m, int zoff, int threads) {
    const int K = g.ncol(), n = g.nrow();
    const int P = seed.size();
    if (charts.size() != K || cf.size() != K || fp.size() != K ||
        u_at.nrow() != n || u_at.ncol() != K) {
        stop("the inputs of the compiled joint filter do not agree in size.");
    }
    const std::size_t nKP = static_cast<std::size_t>(n) * K * P;
    if (U1.size() != static_cast<R_xlen_t>(nKP) ||
        U2.size() != static_cast<R_xlen_t>(nKP * P) ||
        U3.size() != static_cast<R_xlen_t>(nKP * P * P)) {
        stop("'U1', 'U2' and 'U3' must have dimensions c(n, K, P), c(n, K, P, P) and c(n, K, P, P, P).");
    }
    if (nd == 2 && U4.size() != static_cast<R_xlen_t>(nKP * P * P * P)) {
        stop("'U4' must have dimensions c(n, K, P, P, P, P).");
    }
    if (static_cast<int>(dirs.size()) != nd) {
        stop("the number of directions does not agree with the order.");
    }
    for (int d = 0; d < nd; ++d) {
        if (static_cast<int>(dirs[d].size()) != m) {
            stop("a direction must have one value per unknown.");
        }
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
    const double* x0 = dirs[0].data();
    const double* x1 = nd == 2 ? dirs[1].data() : nullptr;

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
        read_cross_pairs(J, f, K, zoff, m);
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
        C.om_t = as<std::vector<double> >(c["om_t"]);
        C.a_t = as<std::vector<double> >(c["a_t"]);
        C.b_t = as<std::vector<double> >(c["b_t"]);
        C.f0_3 = as<std::vector<double> >(c["f0_3"]);
        if (nd == 2) {
            C.om_q = as<std::vector<double> >(c["om_q"]);
            C.a_q = as<std::vector<double> >(c["a_q"]);
            C.b_q = as<std::vector<double> >(c["b_q"]);
            C.f0_4 = as<std::vector<double> >(c["f0_4"]);
        }
        const std::size_t p = J.p, q = J.q, ndz = nd;
        if (C.om_u.size() != static_cast<std::size_t>(m) || C.om_uu.size() != mm ||
            C.f0_u.size() != static_cast<std::size_t>(m) || C.f0_uu.size() != mm ||
            C.a.size() != p || C.a_u.size() != p * m || C.a_uu.size() != p * mm ||
            C.b.size() != q || C.b_u.size() != q * m || C.b_uu.size() != q * mm) {
            stop("the chart of a term does not agree with its orders and the number of unknowns.");
        }
        if (C.om_t.size() != ndz * mm || C.a_t.size() != ndz * p * mm ||
            C.b_t.size() != ndz * q * mm || C.f0_3.size() != ndz * mm) {
            stop("the third derivatives of a term's chart do not agree with its orders and the number of unknowns.");
        }
        if (nd == 2 && (C.om_q.size() != mm || C.a_q.size() != p * mm ||
                        C.b_q.size() != q * mm || C.f0_4.size() != mm)) {
            stop("the fourth derivatives of a term's chart do not agree with its orders and the number of unknowns.");
        }
        // the rings must reach the largest cross lag as well
        L = std::max(L, std::max(std::max(J.p, J.q), J.xmax));
    }
    const int D = L + 1;   // the ring's depth

    // the coefficients of every chart as quantities, with their contractions
    // formed once (the charts are shared by all the groups and threads)
    std::vector<JChartJets> chj(K);
    for (int k = 0; k < K; ++k) {
        const JCurvChart& C = ch[k];
        const JCurvCoef& J = co[k];
        JChartJets& CJ = chj[k];
        const std::size_t p = J.p, q = J.q;
        const std::size_t per = static_cast<std::size_t>(nd) * (m + 1) + 1 + m;
        CJ.store.assign((p + q + 1) * per, 0.0);
        CJ.a.resize(p);
        CJ.b.resize(q);
        auto fill = [&](JetV& V, double* st, double v, const double* d1,
                        const double* d2, const double* t0, const double* t1,
                        const double* d4) {
            V = JetV();
            V.v = v;
            V.d1 = d1;
            V.d2 = d2;
            V.d3[0] = t0;
            V.d3[1] = t1;
            V.d4 = d4;
            double* ad = st;
            double* d2x = st + nd;
            double* d2vw = st + nd + static_cast<std::size_t>(nd) * m;
            double* d3vw = d2vw + 1;
            jet_aux(m, nd, x0, x1, d1, d2, t0, ad, d2x,
                    nd == 2 ? d2vw : nullptr, nd == 2 ? d3vw : nullptr);
            for (int d = 0; d < nd; ++d) {
                V.ad[d] = ad[d];
                V.d2x[d] = d2x + static_cast<std::size_t>(d) * m;
            }
            if (nd == 2) {
                V.d2vw = *d2vw;
                V.d3vw = d3vw;
            }
        };
        for (std::size_t i = 0; i < p; ++i) {
            fill(CJ.a[i], &CJ.store[i * per], C.a[i], &C.a_u[i * m],
                 &C.a_uu[i * mm], &C.a_t[i * mm],
                 nd == 2 ? &C.a_t[(p + i) * mm] : nullptr,
                 nd == 2 ? &C.a_q[i * mm] : nullptr);
        }
        for (std::size_t j = 0; j < q; ++j) {
            fill(CJ.b[j], &CJ.store[(p + j) * per], C.b[j], &C.b_u[j * m],
                 &C.b_uu[j * mm], &C.b_t[j * mm],
                 nd == 2 ? &C.b_t[(q + j) * mm] : nullptr,
                 nd == 2 ? &C.b_q[j * mm] : nullptr);
        }
        fill(CJ.f0, &CJ.store[(p + q) * per], C.f0, C.f0_u.data(),
             C.f0_uu.data(), C.f0_3.data(),
             nd == 2 ? &C.f0_3[mm] : nullptr,
             nd == 2 ? C.f0_4.data() : nullptr);
    }

    NumericMatrix curvature(m, m);
    std::vector<NumericMatrix> jac;
    for (int k = 0; k < K; ++k) jac.push_back(NumericMatrix(n, m));
    std::vector<double*> jac_p(K);
    for (int k = 0; k < K; ++k) jac_p[k] = jac[k].begin();
    // dphi[k][d]: the second derivative of the level applied to direction d;
    // dpsi[k]: the third derivative along the first applied to the second
    std::vector<NumericMatrix> dphi;
    for (int k = 0; k < K; ++k) {
        for (int d = 0; d < nd; ++d) dphi.push_back(NumericMatrix(n, m));
    }
    std::vector<double*> dphi_p(static_cast<std::size_t>(K) * nd);
    for (std::size_t i = 0; i < dphi_p.size(); ++i) dphi_p[i] = dphi[i].begin();
    std::vector<NumericMatrix> dpsi;
    if (nd == 2) {
        for (int k = 0; k < K; ++k) dpsi.push_back(NumericMatrix(n, m));
    }
    std::vector<double*> dpsi_p(dpsi.size());
    for (std::size_t i = 0; i < dpsi.size(); ++i) dpsi_p[i] = dpsi[i].begin();
    const double* gp = g.begin();
    const double* up = u_at.begin();
    const double* U1p = U1.begin();
    const double* U2p = U2.begin();
    const double* U3p = U3.begin();
    const double* U4p = U4.begin();

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
        const std::size_t ns = static_cast<std::size_t>(K) * D;
        const std::size_t Ps = static_cast<std::size_t>(P);
        const std::size_t ms = static_cast<std::size_t>(m);
        JRing F, Uu;
        F.init(ns, m, nd);
        Uu.init(ns, m, nd);
        std::vector<double> ecol(m, 0.0);
        // the predictors at a row: gradients (m x P), and for a filtered
        // equation the Hessian, the third derivatives, the fourth, and the
        // Hessian applied to each direction (null for a static equation)
        std::vector<double> Ed1(Ps * ms);
        std::vector<const double*> Ed2(P), Ed4(P), Ed3(Ps * nd), Ed2x(Ps * nd);
        std::vector<double> TM(ms * Ps), tmp(ms * Ps);
        std::vector<double> dx(Ps * nd), Px(Ps * nd * ms);
        std::vector<double> Psiwv(ms * Ps), phivw(Ps);
        std::vector<double> U1v(Ps), U2m(Ps * Ps), U2x(Ps * nd);
        std::vector<double> U3x(Ps * Ps * nd), XU3(ms * Ps * nd), PU2(ms * Ps);
        std::vector<double> Mm(Ps * Ps), U3wv(Ps), U2phi(Ps);
        // the static starts of the members as quantities: a value, the
        // gradient minus the static row of the member's equation at the
        // group's first row, and every higher derivative zero
        std::vector<double> S1(static_cast<std::size_t>(K) * m, 0.0),
            Z2(mm, 0.0), Zm(ms, 0.0);
        std::vector<JetV> S0(K);
        for (int k = 0; k < K; ++k) {
            if (!co[k].direct || T == 0) continue;
            const int r1 = rows[0] - 1;
            double* s1 = &S1[static_cast<std::size_t>(k) * m];
            for (int i = 0; i < m; ++i) {
                s1[i] = -seed_p[fpe[k]][r1 + static_cast<std::size_t>(n) * i];
            }
            JetV& V = S0[k];
            V.v = co[k].d0[gi];
            V.d1 = s1;
            V.d2 = Z2.data();
            V.d3[0] = Z2.data();
            V.d3[1] = nd == 2 ? Z2.data() : nullptr;
            V.d4 = nd == 2 ? Z2.data() : nullptr;
            const double* xs[2] = {x0, x1};
            for (int d = 0; d < nd; ++d) {
                double s = 0.0;
                for (int i = 0; i < m; ++i) s += s1[i] * xs[d][i];
                V.ad[d] = s;
                V.d2x[d] = Zm.data();
            }
            if (nd == 2) {
                V.d2vw = 0.0;
                V.d3vw = Zm.data();
            }
        }
        for (int t = 0; t < T; ++t) {
            const int ts = t % D;
            // the levels of every filtered parameter at this time
            for (int k = 0; k < K; ++k) {
                const JCurvChart& C = ch[k];
                const JCurvCoef& J = co[k];
                const JChartJets& CJ = chj[k];
                const std::size_t cur = static_cast<std::size_t>(k + K * ts);
                double Jv = C.omega;
                double* J1 = &F.d1[cur * m];
                double* J2 = &F.d2[cur * mm];
                double* J3 = &F.d3[cur * nd * mm];
                double* J4 = nd == 2 ? &F.d4[cur * mm] : nullptr;
                if (J.direct && t == 0) {
                    const double* s1 = &S1[static_cast<std::size_t>(k) * m];
                    std::copy(s1, s1 + m, J1);
                    std::fill(J2, J2 + mm, 0.0);
                    std::fill(J3, J3 + nd * mm, 0.0);
                    if (nd == 2) std::fill(J4, J4 + mm, 0.0);
                    F.v[cur] = J.d0[gi];
                    F.aux(cur, x0, x1);
                    continue;
                }
                std::copy(C.om_u.begin(), C.om_u.end(), J1);
                std::copy(C.om_uu.begin(), C.om_uu.end(), J2);
                std::copy(C.om_t.begin(), C.om_t.end(), J3);
                if (nd == 2) std::copy(C.om_q.begin(), C.om_q.end(), J4);
                for (int i = 1; i <= J.p; ++i) {
                    if (t - i < 0) continue;
                    const int ls = (t - i) % D;
                    // the score loading times the parameter's own lagged score
                    const std::size_t lag = static_cast<std::size_t>(k + K * ls);
                    const JetV A = CJ.a[i - 1];
                    const JetV B = Uu.view(lag);
                    prod_add(m, A.v, A.d1, A.d2, B.v, B.d1, B.d2, &Jv, J1, J2);
                    prod_add3(m, nd, A, B, J3);
                    if (nd == 2) prod_add4(m, A, B, J4);
                }
                // the cross loadings times the partners' lagged scores, one
                // term per pair at its own lag; a cross loading rides the
                // identity link: a one at its own column, no Hessian and no
                // higher derivative
                for (int h = 0; h < J.nx; ++h) {
                    const int xl = J.xlag[h];
                    if (t - xl < 0) continue;
                    const std::size_t lag = static_cast<std::size_t>((J.xpar[h] - 1) + K * ((t - xl) % D));
                    const int col = zoff + J.col_cross[h] - 1;
                    ecol[col] = 1.0;
                    JetV A;
                    A.v = J.cross[h];
                    A.d1 = ecol.data();
                    for (int d = 0; d < nd; ++d) A.ad[d] = dirs[d][col];
                    const JetV B = Uu.view(lag);
                    prod_add(m, A.v, A.d1, nullptr, B.v, B.d1, B.d2, &Jv, J1, J2);
                    prod_add3(m, nd, A, B, J3);
                    if (nd == 2) prod_add4(m, A, B, J4);
                    ecol[col] = 0.0;
                }
                for (int j = 1; j <= J.q; ++j) {
                    const JetV A = CJ.b[j - 1];
                    const JetV B = (t - j >= 0)
                        ? F.view(static_cast<std::size_t>(k + K * ((t - j) % D)))
                        : (J.direct ? S0[k] : CJ.f0);
                    prod_add(m, A.v, A.d1, A.d2, B.v, B.d1, B.d2, &Jv, J1, J2);
                    prod_add3(m, nd, A, B, J3);
                    if (nd == 2) prod_add4(m, A, B, J4);
                }
                F.v[cur] = Jv;
                F.aux(cur, x0, x1);
            }
            const int r = rows[t] - 1;
            // the predictors of every equation at this row
            for (int a = 0; a < P; ++a) {
                for (int i = 0; i < m; ++i) {
                    Ed1[static_cast<std::size_t>(a) * m + i] =
                        seed_p[a][r + static_cast<std::size_t>(n) * i];
                }
                Ed2[a] = nullptr;
                Ed4[a] = nullptr;
                for (int d = 0; d < nd; ++d) {
                    Ed3[static_cast<std::size_t>(a) * nd + d] = nullptr;
                    Ed2x[static_cast<std::size_t>(a) * nd + d] = nullptr;
                }
            }
            for (int k = 0; k < K; ++k) {
                const int a = fpe[k];
                const std::size_t cur = static_cast<std::size_t>(k + K * ts);
                const double* J1 = &F.d1[cur * m];
                double* e1 = &Ed1[static_cast<std::size_t>(a) * m];
                for (int i = 0; i < m; ++i) e1[i] += J1[i];
                Ed2[a] = &F.d2[cur * mm];
                Ed4[a] = nd == 2 ? &F.d4[cur * mm] : nullptr;
                for (int d = 0; d < nd; ++d) {
                    Ed3[static_cast<std::size_t>(a) * nd + d] =
                        &F.d3[(cur * nd + d) * mm];
                    Ed2x[static_cast<std::size_t>(a) * nd + d] =
                        &F.d2x[(cur * nd + d) * m];
                }
                for (int i = 0; i < m; ++i) {
                    jac_p[k][r + static_cast<std::size_t>(n) * i] = e1[i];
                    for (int d = 0; d < nd; ++d) {
                        dphi_p[static_cast<std::size_t>(k) * nd + d]
                              [r + static_cast<std::size_t>(n) * i] =
                            F.d2x[(cur * nd + d) * m + i];
                    }
                    if (nd == 2) {
                        dpsi_p[k][r + static_cast<std::size_t>(n) * i] =
                            F.d3vw[cur * m + i];
                    }
                }
                const double gk = gp[r + static_cast<std::size_t>(n) * k];
                const double* top = nd == 2 ? &F.d4[cur * mm] : &F.d3[cur * nd * mm];
                for (std::size_t i = 0; i < mm; ++i) Wg[i] += gk * top[i];
            }
            // what the driving scores share across parameters: each
            // direction applied to the predictors' gradients and the
            // Hessians applied to each direction (m x P)
            for (int d = 0; d < nd; ++d) {
                for (int a = 0; a < P; ++a) {
                    double s = 0.0;
                    for (int i = 0; i < m; ++i) {
                        s += Ed1[static_cast<std::size_t>(a) * m + i] * dirs[d][i];
                    }
                    dx[static_cast<std::size_t>(d) * P + a] = s;
                    const double* e2x = Ed2x[static_cast<std::size_t>(a) * nd + d];
                    for (int i = 0; i < m; ++i) {
                        Px[(static_cast<std::size_t>(d) * P + a) * m + i] =
                            e2x ? e2x[i] : 0.0;
                    }
                }
            }
            if (nd == 2) {
                // v' Phi_c w, and each predictor's third derivative along w
                // applied to v
                for (int c = 0; c < P; ++c) {
                    double s = 0.0;
                    for (int i = 0; i < m; ++i) {
                        s += Px[static_cast<std::size_t>(c) * m + i] * x1[i];
                    }
                    phivw[c] = s;
                }
                for (int a = 0; a < P; ++a) {
                    const double* e3w = Ed3[static_cast<std::size_t>(a) * nd + 1];
                    for (int i = 0; i < m; ++i) {
                        double s = 0.0;
                        if (e3w) {
                            for (int j = 0; j < m; ++j) {
                                s += e3w[i + static_cast<std::size_t>(m) * j] * x0[j];
                            }
                        }
                        Psiwv[static_cast<std::size_t>(a) * m + i] = s;
                    }
                }
            }
            // the driving score of every filtered parameter at this row
            for (int l = 0; l < K; ++l) {
                const std::size_t cur = static_cast<std::size_t>(l + K * ts);
                const std::size_t s1 = static_cast<std::size_t>(n) * K;
                const std::size_t base = static_cast<std::size_t>(r) + static_cast<std::size_t>(n) * l;
                const double* u1 = U1p + base;
                const double* u2 = U2p + base;
                const double* u3 = U3p + base;
                const double* u4 = nd == 2 ? U4p + base : nullptr;
                Uu.v[cur] = up[base];
                double* d1 = &Uu.d1[cur * m];
                double* d2 = &Uu.d2[cur * mm];
                for (int a = 0; a < P; ++a) U1v[a] = u1[s1 * a];
                for (int b = 0; b < P; ++b) {
                    for (int a = 0; a < P; ++a) {
                        U2m[a + Ps * b] = u2[s1 * (a + Ps * b)];
                    }
                }
                // d1 = Vm U1
                for (int i = 0; i < m; ++i) {
                    double v = 0.0;
                    for (int a = 0; a < P; ++a) {
                        v += Ed1[static_cast<std::size_t>(a) * m + i] * U1v[a];
                    }
                    d1[i] = v;
                }
                // d2 = (Vm U2) Vm' + sum_a U1[a] E_a''
                mat_mulp(m, P, TM.data(), Ed1.data(), U2m.data());
                for (std::size_t i = 0; i < mm; ++i) d2[i] = 0.0;
                add_xyt(m, P, d2, TM.data(), Ed1.data());
                for (int a = 0; a < P; ++a) {
                    if (!Ed2[a]) continue;
                    for (std::size_t i = 0; i < mm; ++i) d2[i] += U1v[a] * Ed2[a][i];
                }
                // d3 along each direction:
                //   Vm U3x Vm' + Px U2 Vm' + Vm U2 Px'
                //   + sum_a (U1[a] E_a''' + (U2 dx)[a] E_a'')
                for (int d = 0; d < nd; ++d) {
                    const double* z = &dx[static_cast<std::size_t>(d) * P];
                    double* X3 = &U3x[Ps * Ps * d];
                    const double* Pxd = &Px[Ps * ms * d];
                    for (int a = 0; a < P; ++a) {
                        double s = 0.0;
                        for (int b = 0; b < P; ++b) s += U2m[a + Ps * b] * z[b];
                        U2x[static_cast<std::size_t>(d) * P + a] = s;
                    }
                    for (int b = 0; b < P; ++b) {
                        for (int a = 0; a < P; ++a) {
                            double s = 0.0;
                            for (int c = 0; c < P; ++c) {
                                s += u3[s1 * (a + Ps * (b + Ps * c))] * z[c];
                            }
                            X3[a + Ps * b] = s;
                        }
                    }
                    mat_mulp(m, P, &XU3[ms * Ps * d], Ed1.data(), X3);
                    mat_mulp(m, P, PU2.data(), Pxd, U2m.data());
                    double* d3o = &Uu.d3[(cur * nd + d) * mm];
                    for (std::size_t i = 0; i < mm; ++i) d3o[i] = 0.0;
                    add_xyt(m, P, d3o, &XU3[ms * Ps * d], Ed1.data());
                    add_xyt(m, P, d3o, PU2.data(), Ed1.data());
                    add_xyt(m, P, d3o, TM.data(), Pxd);
                    for (int a = 0; a < P; ++a) {
                        const double* e3 = Ed3[static_cast<std::size_t>(a) * nd + d];
                        if (e3) {
                            for (std::size_t i = 0; i < mm; ++i) d3o[i] += U1v[a] * e3[i];
                        }
                        if (Ed2[a]) {
                            const double c2 = U2x[static_cast<std::size_t>(d) * P + a];
                            for (std::size_t i = 0; i < mm; ++i) d3o[i] += c2 * Ed2[a][i];
                        }
                    }
                }
                if (nd == 2) {
                    const double* dv = &dx[0];
                    const double* dw = &dx[Ps];
                    const double* Pv = &Px[0];
                    const double* Pw = &Px[Ps * ms];
                    // Mm = U4[., ., dv, dw] + U3[., ., phivw]
                    for (int b = 0; b < P; ++b) {
                        for (int a = 0; a < P; ++a) {
                            double s = 0.0;
                            for (int c = 0; c < P; ++c) {
                                for (int e = 0; e < P; ++e) {
                                    s += u4[s1 * (a + Ps * (b + Ps * (c + Ps * e)))] *
                                         (dv[c] * dw[e]);
                                }
                                s += u3[s1 * (a + Ps * (b + Ps * c))] * phivw[c];
                            }
                            Mm[a + Ps * b] = s;
                        }
                    }
                    for (int a = 0; a < P; ++a) {
                        double s = 0.0, t2 = 0.0;
                        for (int b = 0; b < P; ++b) {
                            s += U3x[Ps * Ps + a + Ps * b] * dv[b];
                            t2 += U2m[a + Ps * b] * phivw[b];
                        }
                        U3wv[a] = s;
                        U2phi[a] = t2;
                    }
                    double* d4o = &Uu.d4[cur * mm];
                    for (std::size_t i = 0; i < mm; ++i) d4o[i] = 0.0;
                    // Vm Mm Vm'
                    mat_mulp(m, P, tmp.data(), Ed1.data(), Mm.data());
                    add_xyt(m, P, d4o, tmp.data(), Ed1.data());
                    // Pw U3x1 Vm' + Vm U3x1 Pw'
                    mat_mulp(m, P, tmp.data(), Pw, &U3x[0]);
                    add_xyt(m, P, d4o, tmp.data(), Ed1.data());
                    add_xyt(m, P, d4o, &XU3[0], Pw);
                    // Pv U3x2 Vm' + Vm U3x2 Pv'
                    mat_mulp(m, P, tmp.data(), Pv, &U3x[Ps * Ps]);
                    add_xyt(m, P, d4o, tmp.data(), Ed1.data());
                    add_xyt(m, P, d4o, &XU3[ms * Ps], Pv);
                    // Psiwv U2 Vm' + Vm U2 Psiwv'
                    mat_mulp(m, P, tmp.data(), Psiwv.data(), U2m.data());
                    add_xyt(m, P, d4o, tmp.data(), Ed1.data());
                    add_xyt(m, P, d4o, TM.data(), Psiwv.data());
                    // Pv U2 Pw' + Pw U2 Pv'
                    mat_mulp(m, P, tmp.data(), Pv, U2m.data());
                    add_xyt(m, P, d4o, tmp.data(), Pw);
                    mat_mulp(m, P, tmp.data(), Pw, U2m.data());
                    add_xyt(m, P, d4o, tmp.data(), Pv);
                    for (int a = 0; a < P; ++a) {
                        if (Ed2[a]) {
                            const double c2 = U3wv[a] + U2phi[a];
                            for (std::size_t i = 0; i < mm; ++i) d4o[i] += c2 * Ed2[a][i];
                        }
                        const double* e3v = Ed3[static_cast<std::size_t>(a) * nd];
                        const double* e3w = Ed3[static_cast<std::size_t>(a) * nd + 1];
                        if (e3w) {
                            const double c3 = U2x[a];
                            for (std::size_t i = 0; i < mm; ++i) d4o[i] += c3 * e3w[i];
                        }
                        if (e3v) {
                            const double c3 = U2x[Ps + a];
                            for (std::size_t i = 0; i < mm; ++i) d4o[i] += c3 * e3v[i];
                        }
                        if (Ed4[a]) {
                            for (std::size_t i = 0; i < mm; ++i) d4o[i] += U1v[a] * Ed4[a][i];
                        }
                    }
                }
                Uu.aux(cur, x0, x1);
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
    List J(K), DP(K);
    for (int k = 0; k < K; ++k) {
        J[k] = jac[k];
        List dk(nd);
        for (int d = 0; d < nd; ++d) dk[d] = dphi[static_cast<std::size_t>(k) * nd + d];
        DP[k] = dk;
    }
    if (nd == 1) {
        return List::create(_["jacobian"] = J, _["curvature"] = curvature,
                            _["dphi"] = DP);
    }
    List PS(K);
    for (int k = 0; k < K; ++k) PS[k] = dpsi[k];
    return List::create(_["jacobian"] = J, _["curvature"] = curvature,
                        _["dphi"] = DP, _["dpsi"] = PS);
}

} // namespace

// The third order: one direction (term_third).
// [[Rcpp::export]]
List gas_joint_third_cpp(List seed, List charts, List cf, IntegerVector fp,
                         NumericMatrix g, NumericMatrix u_at,
                         NumericVector U1, NumericVector U2, NumericVector U3,
                         NumericVector x, List order, int m, int zoff,
                         int threads = 1) {
    std::vector<std::vector<double> > dirs(1);
    dirs[0] = as<std::vector<double> >(x);
    NumericVector none(0);
    return joint_run(1, seed, charts, cf, fp, g, u_at, U1, U2, U3, none, dirs,
                     order, m, zoff, threads);
}

// The fourth order: two directions (term_fourth).
// [[Rcpp::export]]
List gas_joint_fourth_cpp(List seed, List charts, List cf, IntegerVector fp,
                          NumericMatrix g, NumericMatrix u_at,
                          NumericVector U1, NumericVector U2, NumericVector U3,
                          NumericVector U4, NumericVector x, NumericVector y,
                          List order, int m, int zoff, int threads = 1) {
    std::vector<std::vector<double> > dirs(2);
    dirs[0] = as<std::vector<double> >(x);
    dirs[1] = as<std::vector<double> >(y);
    return joint_run(2, seed, charts, cf, fp, g, u_at, U1, U2, U3, U4, dirs,
                     order, m, zoff, threads);
}
