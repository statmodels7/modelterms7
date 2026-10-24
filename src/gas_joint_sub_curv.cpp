#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "gas_par.h"
using namespace Rcpp;

// The second, third and fourth orders of the joint score-driven filter whose
// coefficients vary by observation, compiled: what .gas_joint_sub_core() in
// R/gas_joint_sub_core.R returns, the per-row analogue of
// gas_joint_curv_cpp(), gas_joint_third_cpp() and gas_joint_fourth_cpp() in
// gas_joint_curv.cpp.
//
// The recursion is the one of those kernels: every level and every driving
// score is carried as a value, a gradient, a Hessian, the third derivative
// contracted against each direction and, with two directions, the fourth
// contracted against both, and the products are Leibniz's rule written out
// per order. What changes is the chart. A coefficient of the recursion is no
// longer one set of constants for the whole filter: an omega, a score loading
// or a partial autocorrelation that carries a subformula is
// h(z_t' gamma), read at the row it multiplies. The jets of a coefficient
// are therefore built at every row, on the active columns of the group, from
// the raw design row z, the chained derivative row W = h' z and the link's
// derivatives k2, k3 and k4 (as gas_curvature_gen.cpp reads a value of the
// single-term recursion):
//   d2 = k2 z z',   d3 along x = k3 (z.x) z z',   d4 along x, y = k4 (z.x)(z.y) z z'.
// A scalar coordinate is the same expression with z = 1. The cross loadings
// are values of the recursion like the others (a scalar one has z = 1, a
// developed one its own design row), read at the row they multiply; a value
// whose link has no curvature at any row (the identity link of a cross
// loading) is carried as a value and a gradient only, so a scalar cross
// loading costs what it cost when it was a constant. The autoregressive
// coefficients go through the Levinson-Durbin map, whose derivatives to the
// fourth order are written here as in gas_curvature_gen.cpp. A coefficient
// that does not depend on the row (nothing developed in it) is built once
// per group.
//
// Everything is carried on each group's ACTIVE columns: the unknowns that are
// not parameters of the filter (the fixed effects of the equations), every
// cross loading, and the filter coordinates whose design column is not
// identically zero on the group. With grouping indicators the work per row is
// therefore constant in the number of groups. The curvature of each group is
// stored apart on its active set and added, in group order, to the full
// matrix, so the result does not depend on the thread count.
//
// The starting level of each member is formed at the group's first row, as
// .gas_sub_prep() and .gas_curvature_sub() form it: f0 = omega/(1 - S), with
// S the sum of the autoregressive coefficients, its jets by the expansion of
// the quotient to the third order and by the fixed point
// f0 = omega + S f0 at the fourth.
//
// Every index arriving from R is 1-based. No code in a worker calls the R
// API: a failure is recorded in a flag per group and signalled on the main
// thread.

namespace {

// ---- small dense helpers on m x m column-major blocks -----------------------
inline double dotp(const double* a, const double* b, int m) {
    double s = 0.0;
    for (int i = 0; i < m; ++i) s += a[i] * b[i];
    return s;
}

// out = H v, H symmetric m x m
inline void matvec(const double* H, const double* v, double* out, int m) {
    for (int j = 0; j < m; ++j) out[j] = 0.0;
    for (int k = 0; k < m; ++k) {
        const double vk = v[k];
        for (int j = 0; j < m; ++j) {
            out[j] += H[static_cast<std::size_t>(k) * m + j] * vk;
        }
    }
}

// h += c a b'
inline void add_outer(double* h, double c, const double* a, const double* b,
                      int m) {
    for (int k = 0; k < m; ++k) {
        const double bk = b[k];
        for (int j = 0; j < m; ++j) {
            h[static_cast<std::size_t>(k) * m + j] += c * (a[j] * bk);
        }
    }
}

inline void zero_fill(double* x, std::size_t len) {
    for (std::size_t i = 0; i < len; ++i) x[i] = 0.0;
}

// ---- the inputs --------------------------------------------------------------

// one value of the recursion (the level, a loading, a partial
// autocorrelation) over all the rows: the value, the chained derivative row
// W = h'(eta) z, the raw chart row z (both n x nc, column-major) and the
// link's second to fourth derivatives. A scalar coordinate has nc = 1 and
// z = 1.
struct SVal {
    const double* v = nullptr;
    const double* W = nullptr;
    const double* Z = nullptr;
    const double* k2 = nullptr;
    const double* k3 = nullptr;
    const double* k4 = nullptr;
    int nc = 0;
};

// a member of the joint filter: the orders, the partners and the values of the
// slots omega, kappa_1..p, pacf_1..q and then the cross loadings, the loading
// of partner j (0-based) at lag i (1-based) at the slot 1 + p + q + j + np (i - 1)
struct SMember {
    int p = 0, q = 0;
    std::vector<int> partners;       // 1-based positions among the members
    std::vector<SVal> val;           // 1 + p + q + np p
    std::vector<char> developed;     // the same slots
    std::vector<char> lin;           // the same slots: a cross loading with no
                                     // curvature at any row (value and
                                     // gradient only)
    bool varying_b = false;
};

// where the active coordinates of a value sit: sel the columns of its z, at
// their positions among the group's active unknowns (both 0-based)
struct SGv {
    std::vector<int> at, sel;
};

struct SGroup {
    const int* rows = nullptr;
    int T = 0;
    std::vector<int> act;                      // 0-based unknowns, ascending
    std::vector<std::vector<SGv> > gv;         // [member][slot]
};

struct SInputs {
    int K = 0, n = 0, P = 0, m = 0, zoff = 0, L = 0;
    std::vector<NumericMatrix> seedm;
    std::vector<const double*> seed_p;
    std::vector<int> fpe;
    std::vector<SMember> mem;
    std::vector<SGroup> grp;
    std::vector<IntegerVector> rows_keep;
    std::vector<NumericVector> nv_keep;
    std::vector<NumericMatrix> nm_keep;
};

SVal read_val(List x, int n, int nd, SInputs& I) {
    SVal V;
    NumericVector v = x["v"];
    NumericMatrix W = x["W"];
    NumericMatrix Z = x["Z"];
    NumericVector k2 = x["k2"];
    I.nv_keep.push_back(v);
    I.nm_keep.push_back(W);
    I.nm_keep.push_back(Z);
    I.nv_keep.push_back(k2);
    if (v.size() != n || W.nrow() != n || Z.nrow() != n ||
        W.ncol() != Z.ncol() || k2.size() != n) {
        stop("a value of the joint recursion must have one row per observation.");
    }
    V.v = v.begin();
    V.W = W.begin();
    V.Z = Z.begin();
    V.k2 = k2.begin();
    V.nc = W.ncol();
    if (nd >= 1) {
        NumericVector k3 = x["k3"];
        if (k3.size() != n) stop("'k3' of a value must have one entry per observation.");
        I.nv_keep.push_back(k3);
        V.k3 = k3.begin();
    }
    if (nd >= 2) {
        NumericVector k4 = x["k4"];
        if (k4.size() != n) stop("'k4' of a value must have one entry per observation.");
        I.nv_keep.push_back(k4);
        V.k4 = k4.begin();
    }
    return V;
}

// reads and checks everything but the driving data
void parse_inputs(SInputs& I, List seed, List mem, List groups,
                  IntegerVector fp, int n, int K, int nd, int m, int zoff) {
    I.K = K;
    I.n = n;
    I.m = m;
    I.zoff = zoff;
    I.P = seed.size();
    const int P = I.P;
    if (mem.size() != K || fp.size() != K) {
        stop("the inputs of the compiled joint curvature do not agree in size.");
    }
    I.nv_keep.reserve(static_cast<std::size_t>(K) * 64 + 16);
    I.nm_keep.reserve(static_cast<std::size_t>(K) * 64 + 16);
    I.seedm.resize(P);
    I.seed_p.resize(P);
    for (int a = 0; a < P; ++a) {
        I.seedm[a] = as<NumericMatrix>(seed[a]);
        if (I.seedm[a].nrow() != n || I.seedm[a].ncol() < m) {
            stop("each matrix of 'seed' must have one row per observation and one column per unknown.");
        }
        I.seed_p[a] = I.seedm[a].begin();
    }
    I.fpe.resize(K);
    for (int k = 0; k < K; ++k) {
        I.fpe[k] = fp[k] - 1;
        if (I.fpe[k] < 0 || I.fpe[k] >= P) {
            stop("'fp' must give the equation of each filtered parameter.");
        }
    }
    I.mem.resize(K);
    I.L = 0;
    for (int k = 0; k < K; ++k) {
        List e = mem[k];
        SMember& M = I.mem[k];
        M.p = as<int>(e["p"]);
        M.q = as<int>(e["q"]);
        M.partners = as<std::vector<int> >(e["partners"]);
        const std::size_t np = M.partners.size();
        for (std::size_t j = 0; j < np; ++j) {
            if (M.partners[j] < 1 || M.partners[j] > K) {
                stop("a partner of a filtered parameter is out of range.");
            }
        }
        LogicalVector dev = e["developed"];
        List vl = e["vals"];
        const int nb = 1 + M.p + M.q;
        const int ns = nb + static_cast<int>(np) * M.p;
        if (dev.size() != ns || vl.size() != ns) {
            stop("a term must give one value per base parameter and per cross loading.");
        }
        M.developed.resize(ns);
        for (int s = 0; s < ns; ++s) M.developed[s] = dev[s] ? 1 : 0;
        M.val.resize(ns);
        for (int s = 0; s < ns; ++s) M.val[s] = read_val(as<List>(vl[s]), n, nd, I);
        // a cross loading whose link has no curvature at any row is linear in
        // its coordinates: its second and higher derivatives vanish
        M.lin.assign(ns, 0);
        for (int s = nb; s < ns; ++s) {
            const SVal& V = M.val[s];
            bool zero = true;
            for (int r = 0; r < n && zero; ++r) {
                if (V.k2[r] != 0.0) zero = false;
                if (nd >= 1 && V.k3[r] != 0.0) zero = false;
                if (nd >= 2 && V.k4[r] != 0.0) zero = false;
            }
            M.lin[s] = zero ? 1 : 0;
        }
        M.varying_b = false;
        for (int j = 0; j < M.q; ++j) {
            if (M.developed[1 + M.p + j]) M.varying_b = true;
        }
        I.L = std::max(I.L, std::max(M.p, M.q));
    }
    const int ng = groups.size();
    I.grp.resize(ng);
    I.rows_keep.resize(ng);
    for (int gi = 0; gi < ng; ++gi) {
        List gl = groups[gi];
        SGroup& G = I.grp[gi];
        I.rows_keep[gi] = as<IntegerVector>(gl["rows"]);
        G.rows = I.rows_keep[gi].begin();
        G.T = I.rows_keep[gi].size();
        for (int t = 0; t < G.T; ++t) {
            if (G.rows[t] < 1 || G.rows[t] > n) stop("'order' holds a row outside the data.");
        }
        IntegerVector act = gl["act"];
        const int mk = act.size();
        G.act.resize(mk);
        for (int i = 0; i < mk; ++i) {
            G.act[i] = act[i] - 1;
            if (G.act[i] < 0 || G.act[i] >= m || (i > 0 && G.act[i] <= G.act[i - 1])) {
                stop("the active unknowns of a group must be ascending and within the unknowns.");
            }
        }
        List gv = gl["gv"];
        if (gv.size() != K) stop("a group must describe every term.");
        G.gv.resize(K);
        for (int k = 0; k < K; ++k) {
            List gk = gv[k];
            const SMember& M = I.mem[k];
            const int ns = static_cast<int>(M.val.size());
            if (gk.size() != ns) stop("a group must describe every value of a term.");
            G.gv[k].resize(ns);
            for (int s = 0; s < ns; ++s) {
                List e = gk[s];
                IntegerVector at = e["at"];
                IntegerVector sel = e["sel"];
                if (at.size() != sel.size()) stop("'at' and 'sel' of a value must have one length.");
                SGv& X = G.gv[k][s];
                X.at.resize(at.size());
                X.sel.resize(sel.size());
                for (int i = 0; i < at.size(); ++i) {
                    X.at[i] = at[i] - 1;
                    X.sel[i] = sel[i] - 1;
                    if (X.at[i] < 0 || X.at[i] >= mk || X.sel[i] < 0 ||
                        X.sel[i] >= M.val[s].nc) {
                        stop("a value's active coordinates are out of range.");
                    }
                }
            }
        }
    }
}

// ---- the chart of a value at one row, on a group's active set ---------------

// pointers to the arrays that receive the jets of a coefficient; d3 holds one
// m x m block per direction, and any of d3 and d4 is null at the orders that
// do not need it
struct JOut {
    double* v = nullptr;
    double* d1 = nullptr;
    double* d2 = nullptr;
    double* d3 = nullptr;
    double* d4 = nullptr;
};

// the jets of a value at row r (0-based): h = k2 z z', the third derivative
// along each direction k3 (z.x) z z', the fourth k4 (z.x)(z.y) z z', all on
// the m active unknowns of the group
inline void build_val(const SVal& V, const SGv& G, int r, int n, int m, int nd,
                      const double* x0, const double* x1, double* v,
                      double* d1, double* d2, double* d3, double* d4) {
    const std::size_t mm = static_cast<std::size_t>(m) * m;
    const int na = static_cast<int>(G.at.size());
    const std::size_t rr = static_cast<std::size_t>(r);
    *v = V.v[rr];
    zero_fill(d1, m);
    zero_fill(d2, mm);
    if (d3) zero_fill(d3, static_cast<std::size_t>(nd) * mm);
    if (d4) zero_fill(d4, mm);
    for (int i = 0; i < na; ++i) {
        d1[G.at[i]] = V.W[rr + static_cast<std::size_t>(n) * G.sel[i]];
    }
    const double k2 = V.k2[rr];
    if (k2 != 0.0) {
        for (int j = 0; j < na; ++j) {
            const double zj = V.Z[rr + static_cast<std::size_t>(n) * G.sel[j]];
            if (zj == 0.0) continue;
            for (int i = 0; i < na; ++i) {
                const double zi = V.Z[rr + static_cast<std::size_t>(n) * G.sel[i]];
                d2[G.at[i] + static_cast<std::size_t>(m) * G.at[j]] += k2 * (zi * zj);
            }
        }
    }
    double zd[2] = {0.0, 0.0};
    const double* xs[2] = {x0, x1};
    for (int d = 0; d < nd; ++d) {
        double s = 0.0;
        for (int i = 0; i < na; ++i) {
            s += V.Z[rr + static_cast<std::size_t>(n) * G.sel[i]] * xs[d][G.at[i]];
        }
        zd[d] = s;
    }
    if (d3) {
        for (int d = 0; d < nd; ++d) {
            const double c3 = V.k3[rr] * zd[d];
            if (c3 == 0.0) continue;
            double* T = d3 + static_cast<std::size_t>(d) * mm;
            for (int j = 0; j < na; ++j) {
                const double zj = V.Z[rr + static_cast<std::size_t>(n) * G.sel[j]];
                if (zj == 0.0) continue;
                for (int i = 0; i < na; ++i) {
                    const double zi = V.Z[rr + static_cast<std::size_t>(n) * G.sel[i]];
                    T[G.at[i] + static_cast<std::size_t>(m) * G.at[j]] += c3 * (zi * zj);
                }
            }
        }
    }
    if (d4) {
        const double c4 = V.k4[rr] * zd[0] * zd[1];
        if (c4 != 0.0) {
            for (int j = 0; j < na; ++j) {
                const double zj = V.Z[rr + static_cast<std::size_t>(n) * G.sel[j]];
                if (zj == 0.0) continue;
                for (int i = 0; i < na; ++i) {
                    const double zi = V.Z[rr + static_cast<std::size_t>(n) * G.sel[i]];
                    d4[G.at[i] + static_cast<std::size_t>(m) * G.at[j]] += c4 * (zi * zj);
                }
            }
        }
    }
}

// the value and the gradient of a value with no curvature (all its higher
// jets are zero), at row r on the m active unknowns of the group
inline void build_lin(const SVal& V, const SGv& G, int r, int n, int m,
                      double* v, double* d1) {
    const std::size_t rr = static_cast<std::size_t>(r);
    *v = V.v[rr];
    zero_fill(d1, m);
    const int na = static_cast<int>(G.at.size());
    for (int i = 0; i < na; ++i) {
        d1[G.at[i]] = V.W[rr + static_cast<std::size_t>(n) * G.sel[i]];
    }
}

// ---- the Levinson-Durbin map's derivatives, as gas_levinson2/3/4 ------------
// J is q x q (row j = coefficient, column k = pacf, (j, k) at k q + j); H, T1,
// T2 and Qd hold one q x q block per coefficient, all column-major
struct Lev {
    int q = 0;
    std::vector<double> phi, J, H, T1, T2, Qd;
};

// the second-order map with, optionally, the third along w1 and w2 and the
// fourth along (w1, w2)
inline void levinson(const double* rho, int q, const double* w1,
                     const double* w2, int order, Lev& L) {
    L.q = q;
    const std::size_t qq = static_cast<std::size_t>(q) * q;
    std::vector<double> phi, jac, hes, t1, t2, qd;
    for (int k = 1; k <= q; ++k) {
        std::vector<double> nphi(k), njac(static_cast<std::size_t>(k) * q, 0.0),
            nhes(k * qq, 0.0), nt1(k * qq, 0.0), nt2(k * qq, 0.0),
            nqd(k * qq, 0.0);
        const double rk = rho[k - 1];
        nphi[k - 1] = rk;
        // njac is k x q column-major: (i, c) at c * k + i
        njac[static_cast<std::size_t>(k - 1) * k + (k - 1)] = 1.0;
        const int km1 = k - 1;
        for (int i = 0; i < km1; ++i) {
            const int r = km1 - 1 - i;               // rev index, 0-based
            nphi[i] = phi[i] - rk * phi[r];
            for (int c = 0; c < q; ++c) {
                njac[static_cast<std::size_t>(c) * k + i] =
                    jac[static_cast<std::size_t>(c) * km1 + i] -
                    rk * jac[static_cast<std::size_t>(c) * km1 + r];
            }
            njac[static_cast<std::size_t>(k - 1) * k + i] -= phi[r];
            // the Hessian of coefficient i
            double* h = &nhes[i * qq];
            const double* hi = &hes[i * qq];
            const double* hr = &hes[r * qq];
            for (std::size_t e = 0; e < qq; ++e) h[e] = hi[e] - rk * hr[e];
            for (int c = 0; c < q; ++c) {
                const double jr = jac[static_cast<std::size_t>(c) * km1 + r];
                h[static_cast<std::size_t>(c) * q + (k - 1)] -= jr;     // row k
                h[static_cast<std::size_t>(k - 1) * q + c] -= jr;       // column k
            }
            if (order >= 3) {
                // third along w1 (and w2)
                for (int which = 0; which < ((order >= 4 || w2) ? 2 : 1); ++which) {
                    const double* w = which == 0 ? w1 : w2;
                    if (!w) continue;
                    std::vector<double>& nt = which == 0 ? nt1 : nt2;
                    const std::vector<double>& tt = which == 0 ? t1 : t2;
                    double* T = &nt[i * qq];
                    std::vector<double> hw(q, 0.0);
                    for (int c = 0; c < q; ++c) {
                        for (int b = 0; b < q; ++b) {
                            hw[b] += hr[static_cast<std::size_t>(c) * q + b] * w[c];
                        }
                    }
                    for (std::size_t e = 0; e < qq; ++e) {
                        T[e] = tt[i * qq + e] - rk * tt[r * qq + e] - w[k - 1] * hr[e];
                    }
                    for (int c = 0; c < q; ++c) {
                        T[static_cast<std::size_t>(c) * q + (k - 1)] -= hw[c];
                        T[static_cast<std::size_t>(k - 1) * q + c] -= hw[c];
                    }
                }
            }
            if (order >= 4) {
                double* Qm = &nqd[i * qq];
                const double* t1r = &t1[r * qq];
                const double* t2r = &t2[r * qq];
                std::vector<double> tvw(q, 0.0);
                for (int c = 0; c < q; ++c) {
                    for (int b = 0; b < q; ++b) {
                        tvw[b] += t1r[static_cast<std::size_t>(c) * q + b] * w2[c];
                    }
                }
                for (std::size_t e = 0; e < qq; ++e) {
                    Qm[e] = qd[i * qq + e] - rk * qd[r * qq + e] -
                            w1[k - 1] * t2r[e] - w2[k - 1] * t1r[e];
                }
                for (int c = 0; c < q; ++c) {
                    Qm[static_cast<std::size_t>(c) * q + (k - 1)] -= tvw[c];
                    Qm[static_cast<std::size_t>(k - 1) * q + c] -= tvw[c];
                }
            }
        }
        phi.swap(nphi); jac.swap(njac); hes.swap(nhes);
        t1.swap(nt1); t2.swap(nt2); qd.swap(nqd);
    }
    L.phi = phi;
    L.J = jac;           // q x q, (j, c) at c * q + j
    L.H = hes;
    L.T1 = t1;
    L.T2 = t2;
    L.Qd = qd;
}

// the scratch that the autoregressive coefficients of one member need at one
// row, allocated once per group
struct BScratch {
    int q = 0, m = 0, nd = 0;
    std::vector<double> rho, wdot1, wdot2, hvw, hw, hv4, hw4, hm, t3vw;
    // per partial autocorrelation: its gradient, Hessian, third derivatives
    // (one block per direction), fourth, and the contractions
    std::vector<std::vector<double> > wrows, hrows, trow, qrow, hv1, hv2, tvw;
    Lev Lv, Lm;
    void init(int q_, int m_, int nd_) {
        q = q_;
        m = m_;
        nd = nd_;
        const std::size_t mm = static_cast<std::size_t>(m) * m;
        rho.assign(q, 0.0);
        wdot1.assign(q, 0.0);
        wdot2.assign(q, 0.0);
        hvw.assign(q, 0.0);
        hw.assign(q, 0.0);
        hv4.assign(q, 0.0);
        hw4.assign(q, 0.0);
        hm.assign(q, 0.0);
        t3vw.assign(q, 0.0);
        wrows.assign(q, std::vector<double>(m, 0.0));
        hrows.assign(q, std::vector<double>(mm, 0.0));
        trow.assign(q, std::vector<double>(static_cast<std::size_t>(nd) * mm, 0.0));
        qrow.assign(q, std::vector<double>(nd >= 2 ? mm : 0, 0.0));
        hv1.assign(q, std::vector<double>(m, 0.0));
        hv2.assign(q, std::vector<double>(m, 0.0));
        tvw.assign(q, std::vector<double>(m, 0.0));
    }
};

// The jets of the autoregressive coefficients b_1..b_q of a member at row r:
// b_j = phi_j(rho), rho the vector of the partial autocorrelations, each with
// its own chart. The composition is the one .gas_curvature_sub() writes with
// hb_of(), tb_of() and qb_of(), and gas_curvature_gen.cpp transcribes in
// b_pieces(): the map's own derivatives over the charts' chained rows, plus
// the charts' curvature under the map's Jacobian. out holds the arrays of
// each coefficient (value, gradient, Hessian, third, fourth).
inline void build_b(const SMember& M, const std::vector<SGv>& gvk, int r,
                    int n, int m, int nd, const double* x0, const double* x1,
                    BScratch& S, const JOut* out) {
    const int q = M.q, p = M.p;
    const std::size_t mm = static_cast<std::size_t>(m) * m;
    const std::size_t qq = static_cast<std::size_t>(q) * q;
    for (int k = 0; k < q; ++k) {
        double rk = 0.0;
        build_val(M.val[1 + p + k], gvk[1 + p + k], r, n, m, nd, x0, x1, &rk,
                  S.wrows[k].data(), S.hrows[k].data(),
                  nd >= 1 ? S.trow[k].data() : nullptr,
                  nd >= 2 ? S.qrow[k].data() : nullptr);
        S.rho[k] = rk;
        if (nd >= 1) {
            S.wdot1[k] = dotp(S.wrows[k].data(), x0, m);
            matvec(S.hrows[k].data(), x0, S.hv1[k].data(), m);
        }
        if (nd >= 2) {
            S.wdot2[k] = dotp(S.wrows[k].data(), x1, m);
            matvec(S.hrows[k].data(), x1, S.hv2[k].data(), m);
            S.hvw[k] = dotp(x0, S.hv2[k].data(), m);
            // the third derivative along the first direction, applied to
            // the second
            matvec(S.trow[k].data(), x1, S.tvw[k].data(), m);
        }
    }
    levinson(S.rho.data(), q, nd >= 1 ? S.wdot1.data() : nullptr,
             nd >= 2 ? S.wdot2.data() : nullptr, nd >= 2 ? 4 : (nd == 1 ? 3 : 2),
             S.Lv);
    if (nd >= 2) levinson(S.rho.data(), q, S.hvw.data(), nullptr, 3, S.Lm);
    for (int j = 0; j < q; ++j) {
        const double* Hj = &S.Lv.H[j * qq];
        *out[j].v = S.Lv.phi[j];
        double* d1 = out[j].d1;
        double* h = out[j].d2;
        zero_fill(d1, m);
        zero_fill(h, mm);
        for (int k = 0; k < q; ++k) {
            const double jc = S.Lv.J[static_cast<std::size_t>(k) * q + j];
            if (jc != 0.0) {
                for (int i = 0; i < m; ++i) d1[i] += jc * S.wrows[k][i];
            }
            for (int ll = 0; ll < q; ++ll) {
                const double c = Hj[static_cast<std::size_t>(ll) * q + k];
                if (c != 0.0) add_outer(h, c, S.wrows[k].data(), S.wrows[ll].data(), m);
            }
            if (jc != 0.0) {
                for (std::size_t e = 0; e < mm; ++e) h[e] += jc * S.hrows[k][e];
            }
        }
        for (int d = 0; d < nd; ++d) {
            const std::vector<double>& wd = d == 0 ? S.wdot1 : S.wdot2;
            const std::vector<std::vector<double> >& hvd = d == 0 ? S.hv1 : S.hv2;
            const double* T3 = d == 0 ? &S.Lv.T1[j * qq] : &S.Lv.T2[j * qq];
            double* tbj = out[j].d3 + static_cast<std::size_t>(d) * mm;
            zero_fill(tbj, mm);
            for (int b = 0; b < q; ++b) S.hw[b] = 0.0;
            for (int c = 0; c < q; ++c) {
                for (int b = 0; b < q; ++b) {
                    S.hw[b] += Hj[static_cast<std::size_t>(c) * q + b] * wd[c];
                }
            }
            for (int k = 0; k < q; ++k) {
                for (int ll = 0; ll < q; ++ll) {
                    const double c3 = T3[static_cast<std::size_t>(ll) * q + k];
                    if (c3 != 0.0) add_outer(tbj, c3, S.wrows[k].data(), S.wrows[ll].data(), m);
                    const double c2 = Hj[static_cast<std::size_t>(ll) * q + k];
                    if (c2 != 0.0) {
                        add_outer(tbj, c2, hvd[k].data(), S.wrows[ll].data(), m);
                        add_outer(tbj, c2, S.wrows[k].data(), hvd[ll].data(), m);
                    }
                }
                if (S.hw[k] != 0.0) {
                    for (std::size_t e = 0; e < mm; ++e) tbj[e] += S.hw[k] * S.hrows[k][e];
                }
                const double jc = S.Lv.J[static_cast<std::size_t>(k) * q + j];
                if (jc != 0.0) {
                    const double* tr = S.trow[k].data() + static_cast<std::size_t>(d) * mm;
                    for (std::size_t e = 0; e < mm; ++e) tbj[e] += jc * tr[e];
                }
            }
        }
        if (nd >= 2) {
            // the fourth order of b_j, as qb_of(): the map's fourth derivative
            // and its third along hvw, the charts' pieces in the product rule
            const double* Q4 = &S.Lv.Qd[j * qq];
            const double* T3m = &S.Lm.T1[j * qq];
            const double* T3v = &S.Lv.T1[j * qq];
            const double* T3w = &S.Lv.T2[j * qq];
            double* qbj = out[j].d4;
            zero_fill(qbj, mm);
            for (int b = 0; b < q; ++b) {
                S.hv4[b] = 0.0;
                S.hw4[b] = 0.0;
                S.hm[b] = 0.0;
                S.t3vw[b] = 0.0;
            }
            for (int c = 0; c < q; ++c) {
                for (int b = 0; b < q; ++b) {
                    const double hjbc = Hj[static_cast<std::size_t>(c) * q + b];
                    S.hv4[b] += hjbc * S.wdot1[c];
                    S.hw4[b] += hjbc * S.wdot2[c];
                    S.hm[b] += hjbc * S.hvw[c];
                    S.t3vw[b] += T3v[static_cast<std::size_t>(c) * q + b] * S.wdot2[c];
                }
            }
            for (int k = 0; k < q; ++k) {
                for (int ll = 0; ll < q; ++ll) {
                    const std::size_t kl = static_cast<std::size_t>(ll) * q + k;
                    const double co = Q4[kl] + T3m[kl];
                    if (co != 0.0) add_outer(qbj, co, S.wrows[k].data(), S.wrows[ll].data(), m);
                    if (T3v[kl] != 0.0) {
                        add_outer(qbj, T3v[kl], S.hv2[k].data(), S.wrows[ll].data(), m);
                        add_outer(qbj, T3v[kl], S.wrows[k].data(), S.hv2[ll].data(), m);
                    }
                    if (T3w[kl] != 0.0) {
                        add_outer(qbj, T3w[kl], S.hv1[k].data(), S.wrows[ll].data(), m);
                        add_outer(qbj, T3w[kl], S.wrows[k].data(), S.hv1[ll].data(), m);
                    }
                    const double c2 = Hj[kl];
                    if (c2 != 0.0) {
                        add_outer(qbj, c2, S.tvw[k].data(), S.wrows[ll].data(), m);
                        add_outer(qbj, c2, S.wrows[k].data(), S.tvw[ll].data(), m);
                        add_outer(qbj, c2, S.hv1[k].data(), S.hv2[ll].data(), m);
                        add_outer(qbj, c2, S.hv2[k].data(), S.hv1[ll].data(), m);
                    }
                }
                if (S.t3vw[k] != 0.0) {
                    for (std::size_t e = 0; e < mm; ++e) qbj[e] += S.t3vw[k] * S.hrows[k][e];
                }
                if (S.hm[k] != 0.0) {
                    for (std::size_t e = 0; e < mm; ++e) qbj[e] += S.hm[k] * S.hrows[k][e];
                }
                if (S.hw4[k] != 0.0) {
                    for (std::size_t e = 0; e < mm; ++e) qbj[e] += S.hw4[k] * S.trow[k][e];
                }
                if (S.hv4[k] != 0.0) {
                    const double* tr2 = S.trow[k].data() + mm;
                    for (std::size_t e = 0; e < mm; ++e) qbj[e] += S.hv4[k] * tr2[e];
                }
                const double jc = S.Lv.J[static_cast<std::size_t>(k) * q + j];
                if (jc != 0.0) {
                    for (std::size_t e = 0; e < mm; ++e) qbj[e] += jc * S.qrow[k][e];
                }
            }
        }
    }
}

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

// the merge of the groups' curvature, in group order
inline void scatter_add(NumericMatrix& W, const SGroup& G, const double* Wg,
                        int mfull) {
    const int mk = static_cast<int>(G.act.size());
    double* Wp = W.begin();
    for (int j = 0; j < mk; ++j) {
        for (int i = 0; i < mk; ++i) {
            Wp[G.act[i] + static_cast<std::size_t>(mfull) * G.act[j]] +=
                Wg[i + static_cast<std::size_t>(mk) * j];
        }
    }
}

} // namespace

// ---------------------------------------------------------------------------
// The second order. The coefficients' value, gradient and Hessian are built at
// each row; the recursion is gas_joint_curv_cpp()'s.

// [[Rcpp::export]]
List gas_joint_sub_curv_cpp(List seed, List mem, List groups, IntegerVector fp,
                            NumericMatrix g, NumericMatrix u_at,
                            NumericVector U1, NumericVector U2, int m,
                            int zoff, int threads = 1) {
    const int K = g.ncol(), n = g.nrow();
    if (u_at.nrow() != n || u_at.ncol() != K) {
        stop("the inputs of the compiled joint curvature do not agree in size.");
    }
    SInputs I;
    parse_inputs(I, seed, mem, groups, fp, n, K, 0, m, zoff);
    const int P = I.P;
    if (U1.size() != static_cast<R_xlen_t>(n) * K * P ||
        U2.size() != static_cast<R_xlen_t>(n) * K * P * P) {
        stop("'U1' and 'U2' must have dimensions c(n, K, P) and c(n, K, P, P).");
    }
    const int D = I.L + 1;   // the ring's depth
    const int ng = static_cast<int>(I.grp.size());

    NumericMatrix curvature(m, m);
    std::vector<NumericMatrix> jac;
    for (int k = 0; k < K; ++k) jac.push_back(NumericMatrix(n, m));
    std::vector<double*> jac_p(K);
    for (int k = 0; k < K; ++k) {
        jac_p[k] = jac[k].begin();
        // the columns that are not active on a row hold the static part
        const double* sd = I.seed_p[I.fpe[k]];
        for (std::size_t i = 0; i < static_cast<std::size_t>(n) * m; ++i) jac_p[k][i] = sd[i];
    }
    const double* gp = g.begin();
    const double* up = u_at.begin();
    const double* U1p = U1.begin();
    const double* U2p = U2.begin();
    std::vector<char> bad(ng, 0);

    const int blockg = 256;
    std::vector<std::vector<double> > Wgv;
    std::size_t b0 = 0;

    auto run_group = [&](std::size_t gi) {
        const SGroup& G = I.grp[gi];
        const int* rows = G.rows;
        const int T = G.T;
        const int mk = static_cast<int>(G.act.size());
        const std::size_t mm = static_cast<std::size_t>(mk) * mk;
        double* Wg = Wgv[gi - b0].data();
        // the coefficients of each member: omega, the loadings, the
        // autoregressive coefficients and the cross loadings, as value,
        // gradient and Hessian (no Hessian for a linear cross loading)
        struct C2 {
            double v = 0.0;
            std::vector<double> d1, d2;
        };
        std::vector<std::vector<C2> > cj(K);
        std::vector<C2> f0(K);
        std::vector<BScratch> bs(K);
        for (int k = 0; k < K; ++k) {
            const SMember& M = I.mem[k];
            cj[k].resize(M.val.size());
            for (std::size_t s = 0; s < cj[k].size(); ++s) {
                cj[k][s].d1.assign(mk, 0.0);
                if (!M.lin[s]) cj[k][s].d2.assign(mm, 0.0);
            }
            f0[k].d1.assign(mk, 0.0);
            f0[k].d2.assign(mm, 0.0);
            if (M.q > 0) bs[k].init(M.q, mk, 0);
        }
        auto build_member = [&](int k, int r, bool all) {
            const SMember& M = I.mem[k];
            for (int s = 0; s <= M.p; ++s) {
                if (!all && !M.developed[s]) continue;
                build_val(M.val[s], G.gv[k][s], r, n, mk, 0, nullptr, nullptr,
                          &cj[k][s].v, cj[k][s].d1.data(), cj[k][s].d2.data(),
                          nullptr, nullptr);
            }
            if (M.q > 0 && (all || M.varying_b)) {
                std::vector<JOut> outs(M.q);
                for (int j = 0; j < M.q; ++j) {
                    C2& c = cj[k][1 + M.p + j];
                    outs[j].v = &c.v;
                    outs[j].d1 = c.d1.data();
                    outs[j].d2 = c.d2.data();
                }
                build_b(M, G.gv[k], r, n, mk, 0, nullptr, nullptr, bs[k], outs.data());
            }
            // the cross loadings
            for (int s = 1 + M.p + M.q; s < static_cast<int>(M.val.size()); ++s) {
                if (!all && !M.developed[s]) continue;
                C2& c = cj[k][s];
                if (M.lin[s]) {
                    build_lin(M.val[s], G.gv[k][s], r, n, mk, &c.v, c.d1.data());
                } else {
                    build_val(M.val[s], G.gv[k][s], r, n, mk, 0, nullptr,
                              nullptr, &c.v, c.d1.data(), c.d2.data(), nullptr,
                              nullptr);
                }
            }
        };
        // the starting level of each member from the group's first row
        const int r1 = rows[0] - 1;
        for (int k = 0; k < K; ++k) {
            build_member(k, r1, true);
            const SMember& M = I.mem[k];
            if (M.q == 0) continue;
            double sb = 0.0;
            std::vector<double> Su(mk, 0.0), Suu(mm, 0.0);
            for (int j = 0; j < M.q; ++j) {
                const C2& b = cj[k][1 + M.p + j];
                sb += b.v;
                for (int i = 0; i < mk; ++i) Su[i] += b.d1[i];
                for (std::size_t e = 0; e < mm; ++e) Suu[e] += b.d2[e];
            }
            if (std::fabs(1.0 - sb) < 1e-10) {
                bad[gi] = 1;
                return;
            }
            const double cst = 1.0 / (1.0 - sb);
            const C2& O = cj[k][0];
            f0[k].v = O.v * cst;
            for (int i = 0; i < mk; ++i) {
                f0[k].d1[i] = O.d1[i] * cst + O.v * Su[i] * cst * cst;
            }
            for (int j = 0; j < mk; ++j) {
                for (int i = 0; i < mk; ++i) {
                    const std::size_t ij = i + static_cast<std::size_t>(mk) * j;
                    f0[k].d2[ij] = O.d2[ij] * cst +
                        (O.d1[i] * Su[j] + Su[i] * O.d1[j]) * cst * cst +
                        2.0 * O.v * Su[i] * Su[j] * cst * cst * cst +
                        O.v * Suu[ij] * cst * cst;
                }
            }
        }

        // the rings: level (f) and driving score (u), one slot per
        // (parameter, time mod D)
        const std::size_t ns = static_cast<std::size_t>(K) * D;
        std::vector<double> fv(ns, 0.0), uv(ns, 0.0);
        std::vector<double> fd1(ns * mk, 0.0), ud1(ns * mk, 0.0);
        std::vector<double> fd2(ns * mm, 0.0), ud2(ns * mm, 0.0);
        std::vector<double> Ed1(static_cast<std::size_t>(P) * mk);
        std::vector<const double*> Ed2(P);
        std::vector<double> TM(static_cast<std::size_t>(mk) * P);
        for (int t = 0; t < T; ++t) {
            const int ts = t % D;
            const int r = rows[t] - 1;
            if (t > 0) {
                for (int k = 0; k < K; ++k) build_member(k, r, false);
            }
            // the levels of every filtered parameter at this time
            for (int k = 0; k < K; ++k) {
                const SMember& M = I.mem[k];
                const std::size_t cur = static_cast<std::size_t>(k + K * ts);
                double Jv = cj[k][0].v;
                double* J1 = &fd1[cur * mk];
                double* J2 = &fd2[cur * mm];
                std::copy(cj[k][0].d1.begin(), cj[k][0].d1.end(), J1);
                std::copy(cj[k][0].d2.begin(), cj[k][0].d2.end(), J2);
                const int np = static_cast<int>(M.partners.size());
                for (int i = 1; i <= M.p; ++i) {
                    if (t - i < 0) continue;
                    const int ls = (t - i) % D;
                    // the score loading times the parameter's own lagged score
                    {
                        const std::size_t lag = static_cast<std::size_t>(k + K * ls);
                        const C2& a = cj[k][i];
                        prod_add(mk, a.v, a.d1.data(), a.d2.data(),
                                 uv[lag], &ud1[lag * mk], &ud2[lag * mm],
                                 &Jv, J1, J2);
                    }
                    // the cross loadings times the partners' lagged scores,
                    // each read at this row
                    for (int j = 0; j < np; ++j) {
                        const int l = M.partners[j] - 1;
                        const std::size_t lag = static_cast<std::size_t>(l + K * ls);
                        const std::size_t sx = static_cast<std::size_t>(1 + M.p + M.q) +
                            static_cast<std::size_t>(j) + static_cast<std::size_t>(np) * (i - 1);
                        const C2& a = cj[k][sx];
                        prod_add(mk, a.v, a.d1.data(),
                                 M.lin[sx] ? nullptr : a.d2.data(),
                                 uv[lag], &ud1[lag * mk], &ud2[lag * mm],
                                 &Jv, J1, J2);
                    }
                }
                for (int j = 1; j <= M.q; ++j) {
                    double lv;
                    const double* l1;
                    const double* l2;
                    if (t - j >= 0) {
                        const std::size_t lag = static_cast<std::size_t>(k + K * ((t - j) % D));
                        lv = fv[lag];
                        l1 = &fd1[lag * mk];
                        l2 = &fd2[lag * mm];
                    } else {
                        lv = f0[k].v;
                        l1 = f0[k].d1.data();
                        l2 = f0[k].d2.data();
                    }
                    const C2& b = cj[k][M.p + j];
                    prod_add(mk, b.v, b.d1.data(), b.d2.data(),
                             lv, l1, l2, &Jv, J1, J2);
                }
                fv[cur] = Jv;
            }
            // the predictors of every equation at this row: the static rows,
            // and for a filtered equation the level's gradient on top and
            // its Hessian
            for (int a = 0; a < P; ++a) {
                for (int i = 0; i < mk; ++i) {
                    Ed1[static_cast<std::size_t>(a) * mk + i] =
                        I.seed_p[a][r + static_cast<std::size_t>(n) * G.act[i]];
                }
                Ed2[a] = nullptr;
            }
            for (int k = 0; k < K; ++k) {
                const int a = I.fpe[k];
                const std::size_t cur = static_cast<std::size_t>(k + K * ts);
                const double* J1 = &fd1[cur * mk];
                const double* J2 = &fd2[cur * mm];
                double* e1 = &Ed1[static_cast<std::size_t>(a) * mk];
                for (int i = 0; i < mk; ++i) e1[i] += J1[i];
                Ed2[a] = J2;
                for (int i = 0; i < mk; ++i) {
                    jac_p[k][r + static_cast<std::size_t>(n) * G.act[i]] = e1[i];
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
                double* d1 = &ud1[cur * mk];
                double* d2 = &ud2[cur * mm];
                // d1 = Vm U1
                for (int i = 0; i < mk; ++i) {
                    double v = 0.0;
                    for (int a = 0; a < P; ++a) {
                        v += Ed1[static_cast<std::size_t>(a) * mk + i] * u1[s1 * a];
                    }
                    d1[i] = v;
                }
                // d2 = (Vm U2) Vm' + sum_a U1[a] E_a''
                for (int b = 0; b < P; ++b) {
                    for (int i = 0; i < mk; ++i) {
                        double v = 0.0;
                        for (int a = 0; a < P; ++a) {
                            v += Ed1[static_cast<std::size_t>(a) * mk + i] *
                                 u2[s1 * (a + static_cast<std::size_t>(P) * b)];
                        }
                        TM[i + static_cast<std::size_t>(mk) * b] = v;
                    }
                }
                for (int j = 0; j < mk; ++j) {
                    for (int i = 0; i < mk; ++i) {
                        double v = 0.0;
                        for (int b = 0; b < P; ++b) {
                            v += TM[i + static_cast<std::size_t>(mk) * b] *
                                 Ed1[static_cast<std::size_t>(b) * mk + j];
                        }
                        d2[i + static_cast<std::size_t>(mk) * j] = v;
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
    for (b0 = 0; b0 < static_cast<std::size_t>(ng); b0 += blockg) {
        const std::size_t b1 = std::min(static_cast<std::size_t>(ng),
                                        b0 + static_cast<std::size_t>(blockg));
        Wgv.assign(b1 - b0, std::vector<double>());
        for (std::size_t gi = b0; gi < b1; ++gi) {
            const std::size_t mk = I.grp[gi].act.size();
            Wgv[gi - b0].assign(mk * mk, 0.0);
        }
        if (threads > 1 && ng >= mt7::kMinGroupsPar) {
            RcppParallel::parallelFor(b0, b1, w, 1, threads);
        } else {
            w(b0, b1);
        }
        for (std::size_t gi = b0; gi < b1; ++gi) {
            scatter_add(curvature, I.grp[gi], Wgv[gi - b0].data(), m);
        }
    }
    for (int gi = 0; gi < ng; ++gi) {
        if (bad[gi] == 1) {
            stop("the autoregressive polynomial is at the unit root; the filter has no starting level.");
        }
    }
    // the entries and their transposes collect the same terms in a different
    // order, so the sum is made symmetric exactly
    NumericMatrix sym(m, m);
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            sym(i, j) = 0.5 * (curvature(i, j) + curvature(j, i));
        }
    }
    List J(K);
    for (int k = 0; k < K; ++k) J[k] = jac[k];
    return List::create(_["jacobian"] = J, _["curvature"] = sym);
}

// ---------------------------------------------------------------------------
// The third and fourth orders. Every level and driving score carries, besides
// the value, the gradient and the Hessian, its third derivative contracted
// against each direction (m x m) and, with two directions, its fourth
// derivative contracted against both, as the jets of .gas_joint_core() do.
// The products and the chain rule of the driving scores are those of
// joint_run() in gas_joint_curv.cpp, read on the group's active unknowns, and
// the coefficients' jets are built per row.

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
// a coefficient or a cross loading, whose missing derivatives are null.
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

// the ring of one kind of quantity (levels, driving scores, or the jets of the
// coefficients with one slot): one slot per (parameter, time mod depth), with
// every derivative and the contractions
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
    // a quantity kept as a value and a gradient (and d1 . x_d in ad) only
    void init_lite(int m_, int nd_) {
        m = m_;
        nd = nd_;
        mm = static_cast<std::size_t>(m) * m;
        v.assign(1, 0.0);
        d1.assign(m, 0.0);
        ad.assign(nd, 0.0);
    }
    // the view of a quantity made by init_lite(): every higher derivative null
    JetV lite_view() const {
        JetV V;
        V.v = v[0];
        V.d1 = d1.data();
        for (int d = 0; d < nd; ++d) V.ad[d] = ad[d];
        return V;
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

// the starting level of a member, f0 = omega/(1 - S) with S the sum of its
// autoregressive coefficients, read from the jets of omega (slot 0) and of
// b_1..b_q (slots 1 + p ...) at the group's first row: the third order
// differentiates the quotient, and the fourth is read off the fixed point
// f0 = omega + S f0 (the product rule supplies every term but the one with
// all four derivatives on f0, and dividing by 1 - S removes that one), as
// .gas_chart_prep() forms them. The caller has checked that 1 - S is not
// near zero.
inline void build_f0(const SMember& M, std::vector<JRing>& cjk, int nd,
                     const double* x0, const double* x1, JRing& F0) {
    const int p = M.p, q = M.q, m = F0.m;
    const std::size_t mm = F0.mm;
    const JRing& O = cjk[0];
    const double om = O.v[0];
    const double* om_u = O.d1.data();
    const double* om_uu = O.d2.data();
    double sb = 0.0;
    std::vector<double> Su(m, 0.0), Suu(mm, 0.0);
    std::vector<double> S3(static_cast<std::size_t>(nd) * mm, 0.0);
    std::vector<double> S4(nd == 2 ? mm : 0, 0.0);
    for (int j = 0; j < q; ++j) {
        const JRing& B = cjk[1 + p + j];
        sb += B.v[0];
        for (int i = 0; i < m; ++i) Su[i] += B.d1[i];
        for (std::size_t e = 0; e < mm; ++e) Suu[e] += B.d2[e];
        for (std::size_t e = 0; e < static_cast<std::size_t>(nd) * mm; ++e) S3[e] += B.d3[e];
        if (nd == 2) {
            for (std::size_t e = 0; e < mm; ++e) S4[e] += B.d4[e];
        }
    }
    const double cst = 1.0 / (1.0 - sb);
    F0.v[0] = om * cst;
    double* f1 = F0.d1.data();
    double* f2 = F0.d2.data();
    double* f3 = F0.d3.data();
    for (int i = 0; i < m; ++i) f1[i] = om_u[i] * cst + om * Su[i] * cst * cst;
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            const std::size_t ij = i + static_cast<std::size_t>(m) * j;
            f2[ij] = om_uu[ij] * cst +
                     (om_u[i] * Su[j] + Su[i] * om_u[j]) * cst * cst +
                     2.0 * om * Su[i] * Su[j] * cst * cst * cst +
                     om * Suu[ij] * cst * cst;
        }
    }
    const double* xs[2] = {x0, x1};
    std::vector<double> dS2(m), dom2(m);
    for (int d = 0; d < nd; ++d) {
        const double* vd = xs[d];
        const double dS = dotp(Su.data(), vd, m);
        const double dom = dotp(om_u, vd, m);
        matvec(Suu.data(), vd, dS2.data(), m);
        matvec(om_uu, vd, dom2.data(), m);
        const double* t_om = O.d3.data() + static_cast<std::size_t>(d) * mm;
        const double* S3d = S3.data() + static_cast<std::size_t>(d) * mm;
        double* f3d = f3 + static_cast<std::size_t>(d) * mm;
        const double c2 = cst * cst, c3 = c2 * cst, c4 = c3 * cst;
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < m; ++i) {
                const std::size_t ij = i + static_cast<std::size_t>(m) * j;
                f3d[ij] = t_om[ij] * cst + om_uu[ij] * (c2 * dS) +
                    2.0 * c3 * dS * (om_u[i] * Su[j] + Su[i] * om_u[j]) +
                    c2 * (dom2[i] * Su[j] + om_u[i] * dS2[j] +
                          dS2[i] * om_u[j] + Su[i] * dom2[j]) +
                    (2.0 * dom * c3 + 6.0 * om * c4 * dS) * Su[i] * Su[j] +
                    2.0 * om * c3 * (dS2[i] * Su[j] + Su[i] * dS2[j]) +
                    dom * c2 * Suu[ij] + 2.0 * om * c3 * dS * Suu[ij] +
                    om * c2 * S3d[ij];
            }
        }
    }
    F0.aux(0, x0, x1);
    if (nd == 2) {
        JRing Sr;
        Sr.init(1, m, nd);
        Sr.v[0] = sb;
        std::copy(Su.begin(), Su.end(), Sr.d1.begin());
        std::copy(Suu.begin(), Suu.end(), Sr.d2.begin());
        std::copy(S3.begin(), S3.end(), Sr.d3.begin());
        std::copy(S4.begin(), S4.end(), Sr.d4.begin());
        Sr.aux(0, x0, x1);
        std::vector<double> J4(mm, 0.0);
        // F0.d4 is zero here, as the fixed point requires
        prod_add4(m, Sr.view(0), F0.view(0), J4.data());
        for (std::size_t e = 0; e < mm; ++e) {
            F0.d4[e] = cst * (O.d4[e] + J4[e]);
        }
    }
}

// The recursion for nd = 1 or 2 directions. U3 and U4 are the driving
// scores' third and fourth derivatives in the predictors, arrays of dimension
// c(n, K, P, P, P) and c(n, K, P, P, P, P); U4 is read only with two
// directions.
List joint_sub_run(int nd, List seed, List mem, List groups, IntegerVector fp,
                   NumericMatrix g, NumericMatrix u_at, NumericVector U1,
                   NumericVector U2, NumericVector U3, NumericVector U4,
                   const std::vector<std::vector<double> >& dirs_full,
                   int mfull, int zoff, int threads) {
    const int K = g.ncol(), n = g.nrow();
    if (u_at.nrow() != n || u_at.ncol() != K) {
        stop("the inputs of the compiled joint filter do not agree in size.");
    }
    SInputs I;
    parse_inputs(I, seed, mem, groups, fp, n, K, nd, mfull, zoff);
    const int P = I.P;
    const std::size_t nKP = static_cast<std::size_t>(n) * K * P;
    if (U1.size() != static_cast<R_xlen_t>(nKP) ||
        U2.size() != static_cast<R_xlen_t>(nKP * P) ||
        U3.size() != static_cast<R_xlen_t>(nKP * P * P)) {
        stop("'U1', 'U2' and 'U3' must have dimensions c(n, K, P), c(n, K, P, P) and c(n, K, P, P, P).");
    }
    if (nd == 2 && U4.size() != static_cast<R_xlen_t>(nKP * P * P * P)) {
        stop("'U4' must have dimensions c(n, K, P, P, P, P).");
    }
    if (static_cast<int>(dirs_full.size()) != nd) {
        stop("the number of directions does not agree with the order.");
    }
    for (int d = 0; d < nd; ++d) {
        if (static_cast<int>(dirs_full[d].size()) != mfull) {
            stop("a direction must have one value per unknown.");
        }
    }
    const int D = I.L + 1;   // the ring's depth
    const int ng = static_cast<int>(I.grp.size());

    NumericMatrix curvature(mfull, mfull);
    std::vector<NumericMatrix> jac;
    for (int k = 0; k < K; ++k) jac.push_back(NumericMatrix(n, mfull));
    std::vector<double*> jac_p(K);
    for (int k = 0; k < K; ++k) {
        jac_p[k] = jac[k].begin();
        // the columns that are not active on a row hold the static part
        const double* sd = I.seed_p[I.fpe[k]];
        for (std::size_t i = 0; i < static_cast<std::size_t>(n) * mfull; ++i) jac_p[k][i] = sd[i];
    }
    // dphi[k][d]: the second derivative of the level applied to direction d;
    // dpsi[k]: the third derivative along the first applied to the second
    std::vector<NumericMatrix> dphi;
    for (int k = 0; k < K; ++k) {
        for (int d = 0; d < nd; ++d) dphi.push_back(NumericMatrix(n, mfull));
    }
    std::vector<double*> dphi_p(static_cast<std::size_t>(K) * nd);
    for (std::size_t i = 0; i < dphi_p.size(); ++i) dphi_p[i] = dphi[i].begin();
    std::vector<NumericMatrix> dpsi;
    if (nd == 2) {
        for (int k = 0; k < K; ++k) dpsi.push_back(NumericMatrix(n, mfull));
    }
    std::vector<double*> dpsi_p(dpsi.size());
    for (std::size_t i = 0; i < dpsi.size(); ++i) dpsi_p[i] = dpsi[i].begin();
    const double* gp = g.begin();
    const double* up = u_at.begin();
    const double* U1p = U1.begin();
    const double* U2p = U2.begin();
    const double* U3p = U3.begin();
    const double* U4p = U4.begin();
    std::vector<char> bad(ng, 0);

    const int blockg = 256;
    std::vector<std::vector<double> > Wgv;
    std::size_t b0 = 0;

    auto run_group = [&](std::size_t gi) {
        const SGroup& G = I.grp[gi];
        const int* rows = G.rows;
        const int T = G.T;
        const int m = static_cast<int>(G.act.size());
        const std::size_t mm = static_cast<std::size_t>(m) * m;
        double* Wg = Wgv[gi - b0].data();
        // the directions on the group's active unknowns
        std::vector<std::vector<double> > dirs(nd);
        for (int d = 0; d < nd; ++d) {
            dirs[d].resize(m);
            for (int i = 0; i < m; ++i) dirs[d][i] = dirs_full[d][G.act[i]];
        }
        const double* x0 = dirs[0].data();
        const double* x1 = nd == 2 ? dirs[1].data() : nullptr;
        // the coefficients of each member as quantities with every
        // derivative: omega, the loadings, the autoregressive coefficients and
        // the cross loadings (a linear one keeps its value and gradient only);
        // and the starting level
        std::vector<std::vector<JRing> > cj(K);
        std::vector<JRing> f0r(K);
        std::vector<BScratch> bs(K);
        for (int k = 0; k < K; ++k) {
            const SMember& M = I.mem[k];
            cj[k].resize(M.val.size());
            for (std::size_t s = 0; s < cj[k].size(); ++s) {
                if (M.lin[s]) cj[k][s].init_lite(m, nd);
                else cj[k][s].init(1, m, nd);
            }
            f0r[k].init(1, m, nd);
            if (M.q > 0) bs[k].init(M.q, m, nd);
        }
        auto build_member = [&](int k, int r, bool all) {
            const SMember& M = I.mem[k];
            for (int s = 0; s <= M.p; ++s) {
                if (!all && !M.developed[s]) continue;
                JRing& R = cj[k][s];
                build_val(M.val[s], G.gv[k][s], r, n, m, nd, x0, x1, &R.v[0],
                          R.d1.data(), R.d2.data(), R.d3.data(),
                          nd == 2 ? R.d4.data() : nullptr);
                R.aux(0, x0, x1);
            }
            if (M.q > 0 && (all || M.varying_b)) {
                std::vector<JOut> outs(M.q);
                for (int j = 0; j < M.q; ++j) {
                    JRing& R = cj[k][1 + M.p + j];
                    outs[j].v = &R.v[0];
                    outs[j].d1 = R.d1.data();
                    outs[j].d2 = R.d2.data();
                    outs[j].d3 = R.d3.data();
                    outs[j].d4 = nd == 2 ? R.d4.data() : nullptr;
                }
                build_b(M, G.gv[k], r, n, m, nd, x0, x1, bs[k], outs.data());
                for (int j = 0; j < M.q; ++j) cj[k][1 + M.p + j].aux(0, x0, x1);
            }
            // the cross loadings
            for (int s = 1 + M.p + M.q; s < static_cast<int>(M.val.size()); ++s) {
                if (!all && !M.developed[s]) continue;
                JRing& R = cj[k][s];
                if (M.lin[s]) {
                    build_lin(M.val[s], G.gv[k][s], r, n, m, &R.v[0], R.d1.data());
                    const double* xs[2] = {x0, x1};
                    for (int d = 0; d < nd; ++d) R.ad[d] = dotp(R.d1.data(), xs[d], m);
                } else {
                    build_val(M.val[s], G.gv[k][s], r, n, m, nd, x0, x1, &R.v[0],
                              R.d1.data(), R.d2.data(), R.d3.data(),
                              nd == 2 ? R.d4.data() : nullptr);
                    R.aux(0, x0, x1);
                }
            }
        };
        // the starting level of each member from the group's first row
        const int r1 = rows[0] - 1;
        for (int k = 0; k < K; ++k) {
            build_member(k, r1, true);
            const SMember& M = I.mem[k];
            if (M.q == 0) continue;
            double sb = 0.0;
            for (int j = 0; j < M.q; ++j) sb += cj[k][1 + M.p + j].v[0];
            if (std::fabs(1.0 - sb) < 1e-10) {
                bad[gi] = 1;
                return;
            }
            build_f0(M, cj[k], nd, x0, x1, f0r[k]);
        }

        const std::size_t ns = static_cast<std::size_t>(K) * D;
        const std::size_t Ps = static_cast<std::size_t>(P);
        const std::size_t ms = static_cast<std::size_t>(m);
        JRing F, Uu;
        F.init(ns, m, nd);
        Uu.init(ns, m, nd);
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
        for (int t = 0; t < T; ++t) {
            const int ts = t % D;
            if (t > 0) {
                for (int k = 0; k < K; ++k) build_member(k, rows[t] - 1, false);
            }
            // the levels of every filtered parameter at this time
            for (int k = 0; k < K; ++k) {
                const SMember& M = I.mem[k];
                const std::vector<JRing>& CJ = cj[k];
                const std::size_t cur = static_cast<std::size_t>(k + K * ts);
                double Jv = CJ[0].v[0];
                double* J1 = &F.d1[cur * m];
                double* J2 = &F.d2[cur * mm];
                double* J3 = &F.d3[cur * nd * mm];
                double* J4 = nd == 2 ? &F.d4[cur * mm] : nullptr;
                std::copy(CJ[0].d1.begin(), CJ[0].d1.end(), J1);
                std::copy(CJ[0].d2.begin(), CJ[0].d2.end(), J2);
                std::copy(CJ[0].d3.begin(), CJ[0].d3.end(), J3);
                if (nd == 2) std::copy(CJ[0].d4.begin(), CJ[0].d4.end(), J4);
                const int np = static_cast<int>(M.partners.size());
                for (int i = 1; i <= M.p; ++i) {
                    if (t - i < 0) continue;
                    const int ls = (t - i) % D;
                    // the score loading times the parameter's own lagged score
                    {
                        const std::size_t lag = static_cast<std::size_t>(k + K * ls);
                        const JetV A = CJ[i].view(0);
                        const JetV B = Uu.view(lag);
                        prod_add(m, A.v, A.d1, A.d2, B.v, B.d1, B.d2, &Jv, J1, J2);
                        prod_add3(m, nd, A, B, J3);
                        if (nd == 2) prod_add4(m, A, B, J4);
                    }
                    // the cross loadings times the partners' lagged scores,
                    // each with its jets at this row; one without curvature
                    // has a value and a gradient only
                    for (int j = 0; j < np; ++j) {
                        const int l = M.partners[j] - 1;
                        const std::size_t lag = static_cast<std::size_t>(l + K * ls);
                        const std::size_t sx = static_cast<std::size_t>(1 + M.p + M.q) +
                            static_cast<std::size_t>(j) + static_cast<std::size_t>(np) * (i - 1);
                        const JetV A = M.lin[sx] ? CJ[sx].lite_view() : CJ[sx].view(0);
                        const JetV B = Uu.view(lag);
                        prod_add(m, A.v, A.d1, A.d2, B.v, B.d1, B.d2, &Jv, J1, J2);
                        prod_add3(m, nd, A, B, J3);
                        if (nd == 2) prod_add4(m, A, B, J4);
                    }
                }
                for (int j = 1; j <= M.q; ++j) {
                    const JetV A = CJ[M.p + j].view(0);
                    const JetV B = (t - j >= 0)
                        ? F.view(static_cast<std::size_t>(k + K * ((t - j) % D)))
                        : f0r[k].view(0);
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
                        I.seed_p[a][r + static_cast<std::size_t>(n) * G.act[i]];
                }
                Ed2[a] = nullptr;
                Ed4[a] = nullptr;
                for (int d = 0; d < nd; ++d) {
                    Ed3[static_cast<std::size_t>(a) * nd + d] = nullptr;
                    Ed2x[static_cast<std::size_t>(a) * nd + d] = nullptr;
                }
            }
            for (int k = 0; k < K; ++k) {
                const int a = I.fpe[k];
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
                    const std::size_t at = static_cast<std::size_t>(r) +
                        static_cast<std::size_t>(n) * G.act[i];
                    jac_p[k][at] = e1[i];
                    for (int d = 0; d < nd; ++d) {
                        dphi_p[static_cast<std::size_t>(k) * nd + d][at] =
                            F.d2x[(cur * nd + d) * m + i];
                    }
                    if (nd == 2) dpsi_p[k][at] = F.d3vw[cur * m + i];
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
    for (b0 = 0; b0 < static_cast<std::size_t>(ng); b0 += blockg) {
        const std::size_t b1 = std::min(static_cast<std::size_t>(ng),
                                        b0 + static_cast<std::size_t>(blockg));
        Wgv.assign(b1 - b0, std::vector<double>());
        for (std::size_t gi = b0; gi < b1; ++gi) {
            const std::size_t mk = I.grp[gi].act.size();
            Wgv[gi - b0].assign(mk * mk, 0.0);
        }
        if (threads > 1 && ng >= mt7::kMinGroupsPar) {
            RcppParallel::parallelFor(b0, b1, w, 1, threads);
        } else {
            w(b0, b1);
        }
        for (std::size_t gi = b0; gi < b1; ++gi) {
            scatter_add(curvature, I.grp[gi], Wgv[gi - b0].data(), mfull);
        }
    }
    for (int gi = 0; gi < ng; ++gi) {
        if (bad[gi] == 1) {
            stop("the autoregressive polynomial is at the unit root; the filter has no starting level.");
        }
    }
    // the entries and their transposes collect the same terms in a different
    // order, so the sum is made symmetric exactly
    NumericMatrix sym(mfull, mfull);
    for (int j = 0; j < mfull; ++j) {
        for (int i = 0; i < mfull; ++i) {
            sym(i, j) = 0.5 * (curvature(i, j) + curvature(j, i));
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
        return List::create(_["jacobian"] = J, _["curvature"] = sym,
                            _["dphi"] = DP);
    }
    List PS(K);
    for (int k = 0; k < K; ++k) PS[k] = dpsi[k];
    return List::create(_["jacobian"] = J, _["curvature"] = sym,
                        _["dphi"] = DP, _["dpsi"] = PS);
}

} // namespace

// The third order: one direction (term_third).
// [[Rcpp::export]]
List gas_joint_sub_third_cpp(List seed, List mem, List groups, IntegerVector fp,
                             NumericMatrix g, NumericMatrix u_at,
                             NumericVector U1, NumericVector U2,
                             NumericVector U3, NumericVector x, int m,
                             int zoff, int threads = 1) {
    std::vector<std::vector<double> > dirs(1);
    dirs[0] = as<std::vector<double> >(x);
    NumericVector none(0);
    return joint_sub_run(1, seed, mem, groups, fp, g, u_at, U1, U2, U3, none,
                         dirs, m, zoff, threads);
}

// The fourth order: two directions (term_fourth).
// [[Rcpp::export]]
List gas_joint_sub_fourth_cpp(List seed, List mem, List groups, IntegerVector fp,
                              NumericMatrix g, NumericMatrix u_at,
                              NumericVector U1, NumericVector U2,
                              NumericVector U3, NumericVector U4,
                              NumericVector x, NumericVector y, int m,
                              int zoff, int threads = 1) {
    std::vector<std::vector<double> > dirs(2);
    dirs[0] = as<std::vector<double> >(x);
    dirs[1] = as<std::vector<double> >(y);
    return joint_sub_run(2, seed, mem, groups, fp, g, u_at, U1, U2, U3, U4,
                         dirs, m, zoff, threads);
}
