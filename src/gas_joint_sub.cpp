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
// The context below is the one of gas_joint.cpp, whose copy is private to
// that file.
//
// Every index arriving from R is 1-based.

namespace {

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
// and lam receive the result.
void sym_eigen(int K, double* A, double* Q, double* lam) {
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

struct SubCtx {
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
// covered (see joint_ctx() in gas_joint.cpp for the fields of `fast`).
SubCtx sub_ctx(SEXP fastS, const NumericMatrix& eta_g) {
    SubCtx c;
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
    SubCtx ctx = sub_ctx(fast.isNotNull() ? SEXP(fast.get()) : R_NilValue,
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
