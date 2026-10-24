// The context that the compiled joint score-driven recursions read: the
// family through distributions7's joint registry and the links through
// linkfunctions7's scalar entries, resolved once per call with
// R_GetCCallable, and the symmetric eigen decomposition the scaling uses.
// Shared by gas_joint.cpp (constant coefficients) and gas_joint_sub.cpp
// (coefficients that vary by observation).
#ifndef MODELTERMS7_GAS_JOINT_CTX_H
#define MODELTERMS7_GAS_JOINT_CTX_H

#include <Rcpp.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include <string>
#include <vector>

namespace gasj {

using namespace Rcpp;

typedef int (*lf7_id_t)(const char*);
typedef void (*lf7_inv12p_t)(int, const double*, double, double*, double*,
                             double*);
typedef double (*lf7_clamp_t)(double, double, double);
typedef int (*d7_id_t)(const char*);
typedef int (*d7_cov_t)(int);
typedef int (*d7_safe_t)(int);
typedef void (*d7_joint_t)(int, double, const double*, double*);

constexpr int kJMaxPar = 8;   // the family's parameters and constants
constexpr int kJMaxK = 8;     // the filtered parameters

// The eigen decomposition of a symmetric K x K matrix by cyclic Jacobi
// rotations: A is overwritten, Q (column-major, columns the eigenvectors)
// and lam receive the result. K is small (the filtered parameters), and the
// rotations need no workspace or library.
inline void sym_eigen(int K, double* A, double* Q, double* lam) {
    for (int i = 0; i < K * K; ++i) Q[i] = 0.0;
    for (int i = 0; i < K; ++i) Q[i + K * i] = 1.0;
    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0.0, tot = 0.0;
        for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) {
            const double a = A[i + K * j] * A[i + K * j];
            tot += a;
            if (i != j) off += a;
        }
        if (off <= 1e-32 * tot || off == 0.0) break;
        for (int p = 0; p < K - 1; ++p) for (int q = p + 1; q < K; ++q) {
            const double apq = A[p + K * q];
            if (apq == 0.0) continue;
            const double th = 0.5 * std::atan2(2.0 * apq,
                                               A[q + K * q] - A[p + K * p]);
            const double c = std::cos(th), s = std::sin(th);
            for (int k = 0; k < K; ++k) {
                const double akp = A[k + K * p], akq = A[k + K * q];
                A[k + K * p] = c * akp - s * akq;
                A[k + K * q] = s * akp + c * akq;
            }
            for (int k = 0; k < K; ++k) {
                const double apk = A[p + K * k], aqk = A[q + K * k];
                A[p + K * k] = c * apk - s * aqk;
                A[q + K * k] = s * apk + c * aqk;
            }
            for (int k = 0; k < K; ++k) {
                const double qkp = Q[k + K * p], qkq = Q[k + K * q];
                Q[k + K * p] = c * qkp - s * qkq;
                Q[k + K * q] = s * qkp + c * qkq;
            }
        }
    }
    for (int i = 0; i < K; ++i) lam[i] = A[i + K * i];
}

struct JointCtx {
    bool ok = false;
    bool thread_safe = false;
    int fam = -1, P = 0, nc = 0, K = 0, n = 0;
    int link[kJMaxPar];
    double lpar[kJMaxPar][4];
    double lwr[kJMaxPar], upr[kJMaxPar];
    const double* y = nullptr;
    const double* eta_all = nullptr;   // n x P, the static predictors
    const double* eta_g = nullptr;     // n x K, those of the filtered ones
    const double* cst[kJMaxPar];
    int cst_scalar[kJMaxPar];
    int G[kJMaxK];
    double scaling = 0.0;
    lf7_inv12p_t inv12 = nullptr;
    lf7_clamp_t clamp = nullptr;
    d7_joint_t grad = nullptr, hess = nullptr, info = nullptr,
        dinfo = nullptr;

    // R's `^`, which the R route uses for K = 1
    static inline double rpow(double x, double y) {
        return y == 2.0 ? x * x : R_pow(x, y);
    }

    // the driving vector u (K) and its derivative C (K x K, column-major)
    // at row `row1` and levels f (K); false where the block of the expected
    // information has a non-positive eigenvalue and d = 1/2
    bool step(int row1, const double* f, double* u, double* C) const {
        const int r = row1 - 1;
        double e[kJMaxPar], th[kJMaxPar], h1[kJMaxPar], h2[kJMaxPar];
        for (int a = 0; a < P; ++a) e[a] = eta_all[r + n * a];
        for (int k = 0; k < K; ++k) e[G[k]] = eta_g[r + n * k] + f[k];
        for (int a = 0; a < P; ++a) {
            double h;
            inv12(link[a], lpar[a], e[a], &h, &h1[a], &h2[a]);
            th[a] = clamp(h, lwr[a], upr[a]);
        }
        for (int j = 0; j < nc; ++j) th[P + j] = cst_scalar[j] ? cst[j][0] : cst[j][r];
        const double yi = y[r];
        double g[kJMaxPar], Hp[kJMaxPar * kJMaxPar];
        grad(fam, yi, th, g);
        hess(fam, yi, th, Hp);
        double gl[kJMaxK], HL[kJMaxK * kJMaxK];
        for (int l = 0; l < K; ++l) {
            const int a = G[l];
            gl[l] = g[a] * h1[a];
            for (int k = 0; k < K; ++k) {
                const int b = G[k];
                HL[l + K * k] = Hp[a + P * b] * h1[a] * h1[b];
            }
            HL[l + K * l] += g[a] * h2[a];
        }
        if (scaling == 0.0) {
            for (int l = 0; l < K; ++l) {
                u[l] = gl[l];
                for (int k = 0; k < K; ++k) C[l + K * k] = HL[l + K * k];
            }
            return true;
        }
        double Ep[kJMaxPar * kJMaxPar], dEp[kJMaxPar * kJMaxPar * kJMaxPar];
        info(fam, yi, th, Ep);
        dinfo(fam, yi, th, dEp);
        // the block on the link scale and its derivative in each filtered
        // predictor: I_lk = -E_ab h_a h_b and
        // dI_lk/de_m = -dE_abc h_c h_a h_b - E_ab (d_ac h''_a h_b +
        // d_bc h_a h''_b), with a = G[l], b = G[k], c = G[m]
        double I[kJMaxK * kJMaxK], dI[kJMaxK * kJMaxK * kJMaxK];
        for (int l = 0; l < K; ++l) for (int k = 0; k < K; ++k) {
            const int a = G[l], b = G[k];
            const double Eab = Ep[a + P * b];
            I[l + K * k] = -Eab * h1[a] * h1[b];
            for (int m = 0; m < K; ++m) {
                const int c = G[m];
                double v = -dEp[a + P * b + P * P * c] * h1[c] * h1[a] * h1[b];
                if (l == m) v -= Eab * h2[a] * h1[b];
                if (k == m) v -= Eab * h1[a] * h2[b];
                dI[l + K * k + K * K * m] = v;
            }
        }
        double S[kJMaxK * kJMaxK], dS[kJMaxK * kJMaxK * kJMaxK];
        if (K == 1) {
            const double d = scaling, Iv = I[0];
            S[0] = rpow(Iv, -d);
            dS[0] = -d * rpow(Iv, -d - 1.0) * dI[0];
        } else {
            double A[kJMaxK * kJMaxK], Q[kJMaxK * kJMaxK], lam[kJMaxK],
                rt[kJMaxK];
            for (int i = 0; i < K * K; ++i) A[i] = I[i];
            sym_eigen(K, A, Q, lam);
            for (int i = 0; i < K; ++i) {
                if (scaling == 0.5) {
                    if (!(lam[i] > 0.0)) return false;
                    rt[i] = std::sqrt(lam[i]);
                } else {
                    rt[i] = lam[i];
                }
            }
            for (int a = 0; a < K; ++a) for (int b = 0; b < K; ++b) {
                double v = 0.0;
                for (int i = 0; i < K; ++i) v += Q[a + K * i] * Q[b + K * i] / rt[i];
                S[a + K * b] = v;
            }
            // the derivative of M = I (d = 1) or I^(1/2) (d = 1/2) in the
            // eigenbasis, the latter solving X dX + dX X = dI, and then
            // dS = -S dM S
            for (int m = 0; m < K; ++m) {
                const double* dIm = dI + K * K * m;
                double T[kJMaxK * kJMaxK], W[kJMaxK * kJMaxK],
                    dM[kJMaxK * kJMaxK];
                for (int a = 0; a < K; ++a) for (int j = 0; j < K; ++j) {
                    double v = 0.0;
                    for (int b = 0; b < K; ++b) v += dIm[a + K * b] * Q[b + K * j];
                    W[a + K * j] = v;
                }
                for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) {
                    double v = 0.0;
                    for (int a = 0; a < K; ++a) v += Q[a + K * i] * W[a + K * j];
                    T[i + K * j] = scaling == 0.5 ? v / (rt[i] + rt[j]) : v;
                }
                for (int a = 0; a < K; ++a) for (int j = 0; j < K; ++j) {
                    double v = 0.0;
                    for (int i = 0; i < K; ++i) v += Q[a + K * i] * T[i + K * j];
                    W[a + K * j] = v;
                }
                for (int a = 0; a < K; ++a) for (int b = 0; b < K; ++b) {
                    double v = 0.0;
                    for (int j = 0; j < K; ++j) v += W[a + K * j] * Q[b + K * j];
                    dM[a + K * b] = v;
                }
                for (int a = 0; a < K; ++a) for (int j = 0; j < K; ++j) {
                    double v = 0.0;
                    for (int i = 0; i < K; ++i) v += S[a + K * i] * dM[i + K * j];
                    W[a + K * j] = v;
                }
                for (int a = 0; a < K; ++a) for (int b = 0; b < K; ++b) {
                    double v = 0.0;
                    for (int j = 0; j < K; ++j) v += W[a + K * j] * S[j + K * b];
                    dS[a + K * b + K * K * m] = -v;
                }
            }
        }
        for (int l = 0; l < K; ++l) {
            double v = 0.0;
            for (int c = 0; c < K; ++c) v += S[l + K * c] * gl[c];
            u[l] = v;
            for (int k = 0; k < K; ++k) {
                double w = 0.0;
                for (int c = 0; c < K; ++c) w += S[l + K * c] * HL[c + K * k];
                for (int c = 0; c < K; ++c) w += dS[l + K * c + K * K * k] * gl[c];
                C[l + K * k] = w;
            }
        }
        return true;
    }
};

// The context, or an inert one where a family, a link or a size is not
// covered. `fast` carries: family (the registry's name), constants (a list,
// each of length 1 or n), links (one name per parameter), link_par (a list,
// one numeric vector per parameter), bounds (a list of pairs), y, eta (n x
// P, the static predictors of every parameter), filtered (1-based positions
// of the filtered parameters) and scaling.
inline JointCtx joint_ctx(SEXP fastS, const NumericMatrix& eta_g) {
    JointCtx c;
    if (Rf_isNull(fastS)) return c;
    List fast(fastS);
    d7_id_t d7id = (d7_id_t) R_GetCCallable("distributions7", "d7_scalar_id");
    d7_cov_t cov = (d7_cov_t) R_GetCCallable("distributions7",
                                             "d7_joint_covered");
    lf7_id_t lfid = (lf7_id_t) R_GetCCallable("linkfunctions7",
                                              "lf7_class_id");
    std::string fam = as<std::string>(fast["family"]);
    c.fam = d7id(fam.c_str());
    c.scaling = fast.containsElementNamed("scaling") ?
        as<double>(fast["scaling"]) : 0.0;
    // the registry covers a family for the gradient and the Hessian (2),
    // enough for an unscaled filter, or for all four entries (4)
    if (c.fam < 0 || cov(c.fam) < (c.scaling == 0.0 ? 2 : 4)) return c;
    CharacterVector links = fast["links"];
    List lpar = fast["link_par"];
    List bounds = fast["bounds"];
    List cst = fast["constants"];
    IntegerVector filt = fast["filtered"];
    NumericMatrix eta_all = fast["eta"];
    c.P = links.size();
    c.nc = cst.size();
    c.K = filt.size();
    c.n = eta_all.nrow();
    if (c.P + c.nc > kJMaxPar || c.K > kJMaxK || c.K != eta_g.ncol() ||
        eta_g.nrow() != c.n || eta_all.ncol() != c.P) {
        return c;
    }
    for (int a = 0; a < c.P; ++a) {
        std::string ln = as<std::string>(links[a]);
        c.link[a] = lfid(ln.c_str());
        if (c.link[a] < 0) return c;
        NumericVector lp = lpar[a];
        if (lp.size() > 4) return c;
        for (int j = 0; j < 4; ++j) c.lpar[a][j] = j < lp.size() ? lp[j] : 0.0;
        NumericVector bd = bounds[a];
        c.lwr[a] = bd[0];
        c.upr[a] = bd[1];
    }
    for (int j = 0; j < c.nc; ++j) {
        NumericVector v = cst[j];
        c.cst[j] = v.begin();
        c.cst_scalar[j] = v.size() == 1;
    }
    for (int k = 0; k < c.K; ++k) c.G[k] = filt[k] - 1;
    NumericVector y = fast["y"];
    c.y = y.begin();
    c.eta_all = eta_all.begin();
    c.eta_g = eta_g.begin();
    c.inv12 = (lf7_inv12p_t) R_GetCCallable("linkfunctions7", "lf7_inv12p");
    c.clamp = (lf7_clamp_t) R_GetCCallable("linkfunctions7", "lf7_clamp");
    c.grad = (d7_joint_t) R_GetCCallable("distributions7", "d7_joint_grad");
    c.hess = (d7_joint_t) R_GetCCallable("distributions7", "d7_joint_hess");
    c.info = (d7_joint_t) R_GetCCallable("distributions7", "d7_joint_info");
    c.dinfo = (d7_joint_t) R_GetCCallable("distributions7", "d7_joint_dinfo");
    d7_safe_t safe = (d7_safe_t) R_GetCCallable("distributions7",
                                                "d7_scalar_thread_safe");
    c.thread_safe = safe(c.fam) == 1;
    c.ok = true;
    return c;
}

} // namespace gasj

#endif
