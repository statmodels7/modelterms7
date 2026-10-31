// The fifth order of the score-driven recursion, written out by order: the
// pieces that the compiled kernels gas_fifth_gen.cpp (a term) and
// gas_fifth_joint.cpp (a joint filter) share.
//
// A value of the recursion along the three directions v, w, x carries its
// value, its first derivative d1, its second derivative P2, the third along
// each direction P3[d], the fourth along each pair P4[pr] (pairs vw, vx, wx)
// and the fifth along all three, P5; each matrix has the two free slots.
// The contractions that the product rule reads (a[d] = d1.d, r2[d] = P2 d,
// s2[pr], r3[pr] = P3[d] e, s3, r4 = P4[vw] x) are formed once per value.
//
// The product rule is written out by order: prod2 (4 terms), prod3 (8),
// prod4 (16, as gas_curvature_gen.cpp's prod4) and prod5 (32). A function of
// several inputs (the driving score of the predictors, an autoregressive
// coefficient of the partial autocorrelations) is composed by Faa di Bruno's
// formula, written out by order in compose2 .. compose5 with its named
// contractions of the outer derivatives against the inputs.
#ifndef MT7_GAS_FIFTH_CORE_H
#define MT7_GAS_FIFTH_CORE_H

#include <vector>
#include <cstddef>
#include <algorithm>

namespace mt7f {

typedef std::vector<double> dvec;

// the pairs vw, vx, wx: their two directions and the third
static const int PA[3] = {0, 0, 1};
static const int PB[3] = {1, 2, 2};
static const int PC[3] = {2, 1, 0};
// the pair holding two given directions
inline int pair_of(int d, int e) {
    int lo = d < e ? d : e, hi = d < e ? e : d;
    return lo == 0 ? (hi == 1 ? 0 : 1) : 2;
}

inline void outer_add(double* h, double c, const double* a, const double* b,
                      int mk) {
    if (c == 0.0) return;
    for (int k = 0; k < mk; ++k) {
        double bk = c * b[k];
        if (bk == 0.0) continue;
        for (int j = 0; j < mk; ++j) h[(std::size_t) k * mk + j] += a[j] * bk;
    }
}
// h += a b' + b a'
inline void sym_add(double* h, const double* a, const double* b, int mk) {
    for (int k = 0; k < mk; ++k) {
        double ak = a[k], bk = b[k];
        for (int j = 0; j < mk; ++j) {
            h[(std::size_t) k * mk + j] += a[j] * bk + b[j] * ak;
        }
    }
}
// h += c * X
inline void axpy(double* h, double c, const double* X, std::size_t len) {
    if (c == 0.0) return;
    for (std::size_t e = 0; e < len; ++e) h[e] += c * X[e];
}
inline void mv(const double* H, const double* v, double* out, int mk) {
    for (int j = 0; j < mk; ++j) out[j] = 0.0;
    for (int k = 0; k < mk; ++k) {
        double vk = v[k];
        if (vk == 0.0) continue;
        for (int j = 0; j < mk; ++j) out[j] += H[(std::size_t) k * mk + j] * vk;
    }
}
inline double dt(const double* a, const double* b, int mk) {
    double s = 0.0;
    for (int k = 0; k < mk; ++k) s += a[k] * b[k];
    return s;
}

struct Jet5 {
    int mk = 0;
    bool high = true;     // false: only v and d1 (a predictor of the first order)
    double v = 0.0;
    double a[3] = {0, 0, 0}, s2[3] = {0, 0, 0}, s3 = 0.0;
    dvec buf;
    double *d1, *P2, *P3[3], *P4[3], *P5, *r2[3], *r3[3], *r4;
    void alloc(int mk_) {
        mk = mk_;
        std::size_t m2 = (std::size_t) mk * mk;
        buf.assign(8 * m2 + 8 * (std::size_t) mk, 0.0);
        double* p = buf.data();
        P2 = p; p += m2;
        for (int d = 0; d < 3; ++d) { P3[d] = p; p += m2; }
        for (int d = 0; d < 3; ++d) { P4[d] = p; p += m2; }
        P5 = p; p += m2;
        d1 = p; p += mk;
        for (int d = 0; d < 3; ++d) { r2[d] = p; p += mk; }
        for (int d = 0; d < 3; ++d) { r3[d] = p; p += mk; }
        r4 = p;
    }
    Jet5() {}
    Jet5(const Jet5& o) { alloc(o.mk); copy_from(o); }
    Jet5& operator=(const Jet5& o) {
        if (this != &o) {
            if (mk != o.mk) alloc(o.mk);
            copy_from(o);
        }
        return *this;
    }
    void copy_from(const Jet5& o) {
        high = o.high; v = o.v; s3 = o.s3;
        for (int d = 0; d < 3; ++d) { a[d] = o.a[d]; s2[d] = o.s2[d]; }
        std::copy(o.buf.begin(), o.buf.end(), buf.begin());
    }
    void clear() {
        v = 0.0;
        high = true;
        std::fill(buf.begin(), buf.end(), 0.0);
    }
    // the contractions along the directions
    void finish(const double* const* dirs) {
        for (int d = 0; d < 3; ++d) a[d] = dt(d1, dirs[d], mk);
        if (!high) {
            for (int d = 0; d < 3; ++d) {
                std::fill(r2[d], r2[d] + mk, 0.0);
                std::fill(r3[d], r3[d] + mk, 0.0);
                s2[d] = 0.0;
            }
            std::fill(r4, r4 + mk, 0.0);
            s3 = 0.0;
            return;
        }
        for (int d = 0; d < 3; ++d) mv(P2, dirs[d], r2[d], mk);
        for (int pr = 0; pr < 3; ++pr) {
            s2[pr] = dt(r2[PA[pr]], dirs[PB[pr]], mk);
            mv(P3[PA[pr]], dirs[PB[pr]], r3[pr], mk);
        }
        s3 = dt(r3[0], dirs[2], mk);
        mv(P4[0], dirs[2], r4, mk);
    }
};

// ---- the product rule, one function per order -------------------------------

inline void prod2(const Jet5& A, const Jet5& B, double* out) {
    std::size_t m2 = (std::size_t) A.mk * A.mk;
    axpy(out, A.v, B.P2, m2);
    axpy(out, B.v, A.P2, m2);
    sym_add(out, A.d1, B.d1, A.mk);
}

inline void prod3(const Jet5& A, const Jet5& B, int d, double* out) {
    int mk = A.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    axpy(out, A.v, B.P3[d], m2);
    axpy(out, B.v, A.P3[d], m2);
    axpy(out, A.a[d], B.P2, m2);
    axpy(out, B.a[d], A.P2, m2);
    sym_add(out, A.d1, B.r2[d], mk);
    sym_add(out, A.r2[d], B.d1, mk);
}

inline void prod4(const Jet5& A, const Jet5& B, int pr, double* out) {
    int mk = A.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    int d = PA[pr], e = PB[pr];
    // both directions on B
    axpy(out, A.v, B.P4[pr], m2);
    axpy(out, B.s2[pr], A.P2, m2);
    sym_add(out, A.d1, B.r3[pr], mk);
    // d on A, e on B
    axpy(out, A.a[d], B.P3[e], m2);
    axpy(out, B.a[e], A.P3[d], m2);
    sym_add(out, A.r2[d], B.r2[e], mk);
    // e on A, d on B
    axpy(out, A.a[e], B.P3[d], m2);
    axpy(out, B.a[d], A.P3[e], m2);
    sym_add(out, A.r2[e], B.r2[d], mk);
    // both on A
    axpy(out, A.s2[pr], B.P2, m2);
    axpy(out, B.v, A.P4[pr], m2);
    sym_add(out, A.r3[pr], B.d1, mk);
}

inline void prod5(const Jet5& A, const Jet5& B, double* out) {
    int mk = A.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    // v, w, x on B
    axpy(out, A.v, B.P5, m2);
    axpy(out, B.s3, A.P2, m2);
    sym_add(out, A.d1, B.r4, mk);
    // one direction on A, the pair of the other two on B
    for (int d = 0; d < 3; ++d) {
        int pr = d == 0 ? 2 : (d == 1 ? 1 : 0);   // the pair without d
        axpy(out, A.a[d], B.P4[pr], m2);
        axpy(out, B.s2[pr], A.P3[d], m2);
        sym_add(out, A.r2[d], B.r3[pr], mk);
    }
    // a pair on A, the third direction on B
    for (int pr = 0; pr < 3; ++pr) {
        int f = PC[pr];
        axpy(out, A.s2[pr], B.P3[f], m2);
        axpy(out, B.a[f], A.P4[pr], m2);
        sym_add(out, A.r3[pr], B.r2[f], mk);
    }
    // v, w, x on A
    axpy(out, A.s3, B.P2, m2);
    axpy(out, B.v, A.P5, m2);
    sym_add(out, A.r4, B.d1, mk);
}

// C += A B at every order (the contractions of C are left to finish())
inline void prod_into(const Jet5& A, const Jet5& B, Jet5& C) {
    C.v += A.v * B.v;
    for (int k = 0; k < A.mk; ++k) C.d1[k] += A.v * B.d1[k] + B.v * A.d1[k];
    prod2(A, B, C.P2);
    for (int d = 0; d < 3; ++d) prod3(A, B, d, C.P3[d]);
    for (int pr = 0; pr < 3; ++pr) prod4(A, B, pr, C.P4[pr]);
    prod5(A, B, C.P5);
}

// the jet of h(z'u): value k0, first derivative k1 z, and k_n times the
// directions' components along z times z z'
inline void rank_one(const double* K, const double* z, const double* zd,
                     Jet5& J) {
    int mk = J.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    J.high = true;
    J.v = K[0];
    for (int k = 0; k < mk; ++k) J.d1[k] = K[1] * z[k];
    std::fill(J.P2, J.P2 + m2, 0.0);
    outer_add(J.P2, K[2], z, z, mk);
    for (int d = 0; d < 3; ++d) {
        std::fill(J.P3[d], J.P3[d] + m2, 0.0);
        outer_add(J.P3[d], K[3] * zd[d], z, z, mk);
    }
    for (int pr = 0; pr < 3; ++pr) {
        std::fill(J.P4[pr], J.P4[pr] + m2, 0.0);
        outer_add(J.P4[pr], K[4] * zd[PA[pr]] * zd[PB[pr]], z, z, mk);
    }
    std::fill(J.P5, J.P5 + m2, 0.0);
    outer_add(J.P5, K[5] * zd[0] * zd[1] * zd[2], z, z, mk);
}

// the starting level, read off the fixed point f0 = omega + S f0: at each
// order (1 - S)^(-1) times the level's derivative plus the product rule of
// S f0 without the term that carries every derivative on f0, which is still
// zero when the order is formed
inline void start_level(const Jet5& om, const Jet5& S, Jet5& f0,
                        const double* const* dirs) {
    int mk = om.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    double c = 1.0 / (1.0 - S.v);
    f0.clear();
    f0.v = om.v * c;
    for (int k = 0; k < mk; ++k) f0.d1[k] = c * (om.d1[k] + S.d1[k] * f0.v);
    f0.finish(dirs);
    dvec acc(m2);
    std::fill(acc.begin(), acc.end(), 0.0);
    prod2(S, f0, acc.data());
    for (std::size_t e = 0; e < m2; ++e) f0.P2[e] = c * (om.P2[e] + acc[e]);
    f0.finish(dirs);
    for (int d = 0; d < 3; ++d) {
        std::fill(acc.begin(), acc.end(), 0.0);
        prod3(S, f0, d, acc.data());
        for (std::size_t e = 0; e < m2; ++e) f0.P3[d][e] = c * (om.P3[d][e] + acc[e]);
    }
    f0.finish(dirs);
    for (int pr = 0; pr < 3; ++pr) {
        std::fill(acc.begin(), acc.end(), 0.0);
        prod4(S, f0, pr, acc.data());
        for (std::size_t e = 0; e < m2; ++e) f0.P4[pr][e] = c * (om.P4[pr][e] + acc[e]);
    }
    f0.finish(dirs);
    std::fill(acc.begin(), acc.end(), 0.0);
    prod5(S, f0, acc.data());
    for (std::size_t e = 0; e < m2; ++e) f0.P5[e] = c * (om.P5[e] + acc[e]);
    f0.finish(dirs);
}

// ---- the composition, one function per order --------------------------------
//
// The function u(e_0, .., e_{P-1}) of P inputs, with outer derivatives
// u_{r_1 .. r_k} given by `get(k, idx)`, k = 1 .. 5. The contraction
// u[.; y_1 .. y_j] of order k against j vectors over the inputs leaves k - j
// indices free; `ctr` forms it with 0, 1 or 2 free indices (the free ones
// first). The vectors over the inputs: g[d] = (V_r . d), y2[pr] = the inputs'
// second derivatives along the pair (zero for a first-order input), y3 =
// their third along v, w, x.

template <typename Get>
struct Outer {
    const Get* get;
    int P;
    // sum over the contracted indices of u_{free.., fixed..} prod y_i
    void ctr(int k, int nfree, int nfix, const double* const* ys, double* out) const {
        int nout = nfree == 0 ? 1 : (nfree == 1 ? P : P * P);
        for (int o = 0; o < nout; ++o) out[o] = 0.0;
        int idx[5];
        int total = 1;
        for (int i = 0; i < k; ++i) total *= P;
        for (int c = 0; c < total; ++c) {
            int rem = c;
            for (int i = k - 1; i >= 0; --i) { idx[i] = rem % P; rem /= P; }
            double w = 1.0;
            for (int i = 0; i < nfix; ++i) {
                w *= ys[i][idx[nfree + i]];
                if (w == 0.0) break;
            }
            if (w == 0.0) continue;
            int o = nfree == 0 ? 0 : (nfree == 1 ? idx[0] : idx[0] + P * idx[1]);
            out[o] += w * (*get)(k, idx);
        }
    }
};

// the workspace of a composition: the vectors over the inputs and the
// coefficient arrays, allocated once
struct CompWork {
    int P = 0, mk = 0;
    dvec g, y2, y3, c1, cv, cm, tmpv, tmpm, row, acc;
    void alloc(int P_, int mk_) {
        P = P_; mk = mk_;
        g.assign(3 * (std::size_t) P, 0.0);
        y2.assign(3 * (std::size_t) P, 0.0);
        y3.assign(P, 0.0);
        c1.assign(P, 0.0);
        cv.assign(P, 0.0);
        cm.assign((std::size_t) P * P, 0.0);
        tmpv.assign(P, 0.0);
        tmpm.assign((std::size_t) P * P, 0.0);
        row.assign(mk, 0.0);
        acc.assign((std::size_t) P * mk, 0.0);
    }
};

// out += V' C V, C a P x P coefficient matrix over the inputs' first
// derivatives
inline void add_VCV(const std::vector<const Jet5*>& E, const double* C, int P,
                    int mk, double* out, dvec& row) {
    for (int rb = 0; rb < P; ++rb) {
        std::fill(row.begin(), row.end(), 0.0);
        bool any = false;
        for (int ra = 0; ra < P; ++ra) {
            double c = C[ra + (std::size_t) P * rb];
            if (c == 0.0) continue;
            const double* Va = E[ra]->d1;
            for (int k = 0; k < mk; ++k) row[k] += c * Va[k];
            any = true;
        }
        if (any) outer_add(out, 1.0, row.data(), E[rb]->d1, mk);
    }
}

// sum_rb sym( sum_ra C[ra, rb] X_ra, Y_rb ), with X_ra a row of input ra (its
// first derivative, or a row of a higher-order input) and Y_rb a row of a
// higher-order input
template <typename RowX, typename RowY>
inline void add_sym_rows(const std::vector<const Jet5*>& E, const double* C,
                         int P, int mk, RowX rx, RowY ry, double* out,
                         dvec& row) {
    for (int rb = 0; rb < P; ++rb) {
        if (!E[rb]->high) continue;
        std::fill(row.begin(), row.end(), 0.0);
        bool any = false;
        for (int ra = 0; ra < P; ++ra) {
            double c = C[ra + (std::size_t) P * rb];
            if (c == 0.0) continue;
            const double* X = rx(E[ra]);
            if (!X) continue;
            for (int k = 0; k < mk; ++k) row[k] += c * X[k];
            any = true;
        }
        if (any) sym_add(out, row.data(), ry(E[rb]), mk);
    }
}

template <typename Get>
inline void comp_vectors(const std::vector<const Jet5*>& E, CompWork& W) {
    int P = W.P;
    for (int r = 0; r < P; ++r) {
        for (int d = 0; d < 3; ++d) {
            W.g[(std::size_t) d * P + r] = E[r]->a[d];
            W.y2[(std::size_t) d * P + r] = E[r]->high ? E[r]->s2[d] : 0.0;
        }
        W.y3[r] = E[r]->high ? E[r]->s3 : 0.0;
    }
}

// order 2: sum_r u_r Phi_r + V' u.. V
template <typename Get>
inline void compose2(const Outer<Get>& U, const std::vector<const Jet5*>& E,
                     CompWork& W, double* out) {
    int P = W.P, mk = W.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    std::fill(out, out + m2, 0.0);
    U.ctr(1, 1, 0, nullptr, W.c1.data());
    for (int r = 0; r < P; ++r) if (E[r]->high) axpy(out, W.c1[r], E[r]->P2, m2);
    U.ctr(2, 2, 0, nullptr, W.cm.data());
    add_VCV(E, W.cm.data(), P, mk, out, W.row);
}

// order 3 along d: the inputs' third derivatives, their second times
// u_{r, g_d}, V' u_{.., g_d} V, and the first derivatives against the second
template <typename Get>
inline void compose3(const Outer<Get>& U, const std::vector<const Jet5*>& E,
                     CompWork& W, int d, double* out) {
    int P = W.P, mk = W.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    std::fill(out, out + m2, 0.0);
    const double* gd = &W.g[(std::size_t) d * P];
    for (int r = 0; r < P; ++r) if (E[r]->high) axpy(out, W.c1[r], E[r]->P3[d], m2);
    const double* ys[3] = {gd};
    U.ctr(2, 1, 1, ys, W.cv.data());
    for (int r = 0; r < P; ++r) if (E[r]->high) axpy(out, W.cv[r], E[r]->P2, m2);
    U.ctr(3, 2, 1, ys, W.tmpm.data());
    add_VCV(E, W.tmpm.data(), P, mk, out, W.row);
    U.ctr(2, 2, 0, nullptr, W.cm.data());
    add_sym_rows(E, W.cm.data(), P, mk,
                 [](const Jet5* e) { return (const double*) e->d1; },
                 [d](const Jet5* e) { return (const double*) e->r2[d]; },
                 out, W.row);
}

// order 4 along the pair pr = (d, e)
template <typename Get>
inline void compose4(const Outer<Get>& U, const std::vector<const Jet5*>& E,
                     CompWork& W, int pr, double* out) {
    int P = W.P, mk = W.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    int d = PA[pr], e = PB[pr];
    std::fill(out, out + m2, 0.0);
    const double* gd = &W.g[(std::size_t) d * P];
    const double* ge = &W.g[(std::size_t) e * P];
    const double* yde = &W.y2[(std::size_t) pr * P];
    // the free slots in one block on input r
    for (int r = 0; r < P; ++r) if (E[r]->high) axpy(out, W.c1[r], E[r]->P4[pr], m2);
    {
        const double* ys[1] = {ge};
        U.ctr(2, 1, 1, ys, W.cv.data());
        for (int r = 0; r < P; ++r) if (E[r]->high) axpy(out, W.cv[r], E[r]->P3[d], m2);
        ys[0] = gd;
        U.ctr(2, 1, 1, ys, W.cv.data());
        for (int r = 0; r < P; ++r) if (E[r]->high) axpy(out, W.cv[r], E[r]->P3[e], m2);
    }
    {
        const double* ys2[2] = {gd, ge};
        U.ctr(3, 1, 2, ys2, W.cv.data());
        const double* ys1[1] = {yde};
        U.ctr(2, 1, 1, ys1, W.tmpv.data());
        for (int r = 0; r < P; ++r) {
            if (E[r]->high) axpy(out, W.cv[r] + W.tmpv[r], E[r]->P2, m2);
        }
    }
    // the free slots in two blocks of first derivatives
    {
        const double* ys2[2] = {gd, ge};
        U.ctr(4, 2, 2, ys2, W.cm.data());
        const double* ys1[1] = {yde};
        U.ctr(3, 2, 1, ys1, W.tmpm.data());
        for (std::size_t o = 0; o < W.cm.size(); ++o) W.cm[o] += W.tmpm[o];
        add_VCV(E, W.cm.data(), P, mk, out, W.row);
    }
    // a first derivative against a second along d, and along e
    {
        const double* ys1[1] = {ge};
        U.ctr(3, 2, 1, ys1, W.cm.data());
        add_sym_rows(E, W.cm.data(), P, mk,
                     [](const Jet5* x) { return (const double*) x->d1; },
                     [d](const Jet5* x) { return (const double*) x->r2[d]; },
                     out, W.row);
        ys1[0] = gd;
        U.ctr(3, 2, 1, ys1, W.cm.data());
        add_sym_rows(E, W.cm.data(), P, mk,
                     [](const Jet5* x) { return (const double*) x->d1; },
                     [e](const Jet5* x) { return (const double*) x->r2[e]; },
                     out, W.row);
    }
    // a first derivative against a third along the pair, and a second along
    // d against a second along e
    U.ctr(2, 2, 0, nullptr, W.cm.data());
    add_sym_rows(E, W.cm.data(), P, mk,
                 [](const Jet5* x) { return (const double*) x->d1; },
                 [pr](const Jet5* x) { return (const double*) x->r3[pr]; },
                 out, W.row);
    add_sym_rows(E, W.cm.data(), P, mk,
                 [d](const Jet5* x) { return x->high ? (const double*) x->r2[d] : (const double*) nullptr; },
                 [e](const Jet5* x) { return (const double*) x->r2[e]; },
                 out, W.row);
}

// order 5 along v, w, x
template <typename Get>
inline void compose5(const Outer<Get>& U, const std::vector<const Jet5*>& E,
                     CompWork& W, double* out) {
    int P = W.P, mk = W.mk;
    std::size_t m2 = (std::size_t) mk * mk;
    std::fill(out, out + m2, 0.0);
    const double* g[3] = {&W.g[0], &W.g[(std::size_t) P], &W.g[(std::size_t) 2 * P]};
    const double* y2[3] = {&W.y2[0], &W.y2[(std::size_t) P], &W.y2[(std::size_t) 2 * P]};
    const double* y3 = W.y3.data();
    // case one: the free slots in one block on input r
    for (int r = 0; r < P; ++r) if (E[r]->high) axpy(out, W.c1[r], E[r]->P5, m2);
    for (int pr = 0; pr < 3; ++pr) {
        const double* ys[1] = {g[PC[pr]]};
        U.ctr(2, 1, 1, ys, W.cv.data());
        for (int r = 0; r < P; ++r) if (E[r]->high) axpy(out, W.cv[r], E[r]->P4[pr], m2);
    }
    for (int d = 0; d < 3; ++d) {
        int pr = d == 0 ? 2 : (d == 1 ? 1 : 0);     // the other two
        const double* ys2[2] = {g[PA[pr]], g[PB[pr]]};
        U.ctr(3, 1, 2, ys2, W.cv.data());
        const double* ys1[1] = {y2[pr]};
        U.ctr(2, 1, 1, ys1, W.tmpv.data());
        for (int r = 0; r < P; ++r) {
            if (E[r]->high) axpy(out, W.cv[r] + W.tmpv[r], E[r]->P3[d], m2);
        }
    }
    {
        const double* ys3[3] = {g[0], g[1], g[2]};
        U.ctr(4, 1, 3, ys3, W.cv.data());
        for (int pr = 0; pr < 3; ++pr) {
            const double* ys2[2] = {y2[pr], g[PC[pr]]};
            U.ctr(3, 1, 2, ys2, W.tmpv.data());
            for (int r = 0; r < P; ++r) W.cv[r] += W.tmpv[r];
        }
        const double* ys1[1] = {y3};
        U.ctr(2, 1, 1, ys1, W.tmpv.data());
        for (int r = 0; r < P; ++r) {
            if (E[r]->high) axpy(out, W.cv[r] + W.tmpv[r], E[r]->P2, m2);
        }
    }
    // case two (i): two blocks of first derivatives
    {
        const double* ys3[3] = {g[0], g[1], g[2]};
        U.ctr(5, 2, 3, ys3, W.cm.data());
        for (int pr = 0; pr < 3; ++pr) {
            const double* ys2[2] = {y2[pr], g[PC[pr]]};
            U.ctr(4, 2, 2, ys2, W.tmpm.data());
            for (std::size_t o = 0; o < W.cm.size(); ++o) W.cm[o] += W.tmpm[o];
        }
        const double* ys1[1] = {y3};
        U.ctr(3, 2, 1, ys1, W.tmpm.data());
        for (std::size_t o = 0; o < W.cm.size(); ++o) W.cm[o] += W.tmpm[o];
        add_VCV(E, W.cm.data(), P, mk, out, W.row);
    }
    // (ii): a first derivative against a second along d, the other two
    // directions contracted
    for (int d = 0; d < 3; ++d) {
        int pr = d == 0 ? 2 : (d == 1 ? 1 : 0);
        const double* ys2[2] = {g[PA[pr]], g[PB[pr]]};
        U.ctr(4, 2, 2, ys2, W.cm.data());
        const double* ys1[1] = {y2[pr]};
        U.ctr(3, 2, 1, ys1, W.tmpm.data());
        for (std::size_t o = 0; o < W.cm.size(); ++o) W.cm[o] += W.tmpm[o];
        add_sym_rows(E, W.cm.data(), P, mk,
                     [](const Jet5* x) { return (const double*) x->d1; },
                     [d](const Jet5* x) { return (const double*) x->r2[d]; },
                     out, W.row);
    }
    // (iii): a first derivative against a third along a pair, the third
    // direction contracted
    for (int pr = 0; pr < 3; ++pr) {
        const double* ys1[1] = {g[PC[pr]]};
        U.ctr(3, 2, 1, ys1, W.cm.data());
        add_sym_rows(E, W.cm.data(), P, mk,
                     [](const Jet5* x) { return (const double*) x->d1; },
                     [pr](const Jet5* x) { return (const double*) x->r3[pr]; },
                     out, W.row);
    }
    // (iv): a first derivative against a fourth along v, w, x
    U.ctr(2, 2, 0, nullptr, W.cm.data());
    add_sym_rows(E, W.cm.data(), P, mk,
                 [](const Jet5* x) { return (const double*) x->d1; },
                 [](const Jet5* x) { return (const double*) x->r4; },
                 out, W.row);
    // (vi): a second along d against a third along the other pair
    for (int d = 0; d < 3; ++d) {
        int pr = d == 0 ? 2 : (d == 1 ? 1 : 0);
        add_sym_rows(E, W.cm.data(), P, mk,
                     [d](const Jet5* x) { return x->high ? (const double*) x->r2[d] : (const double*) nullptr; },
                     [pr](const Jet5* x) { return (const double*) x->r3[pr]; },
                     out, W.row);
    }
    // (v): a second along d against a second along e, the third direction
    // contracted
    for (int pr = 0; pr < 3; ++pr) {
        int d = PA[pr], e = PB[pr];
        const double* ys1[1] = {g[PC[pr]]};
        U.ctr(3, 2, 1, ys1, W.cm.data());
        add_sym_rows(E, W.cm.data(), P, mk,
                     [d](const Jet5* x) { return x->high ? (const double*) x->r2[d] : (const double*) nullptr; },
                     [e](const Jet5* x) { return (const double*) x->r2[e]; },
                     out, W.row);
    }
}

// the composed jet: value given, first derivative sum_r u_r V_r, and every
// order
template <typename Get>
inline void compose(const Outer<Get>& U, const std::vector<const Jet5*>& E,
                    double value, CompWork& W, Jet5& out,
                    const double* const* dirs) {
    int P = W.P, mk = W.mk;
    comp_vectors<Get>(E, W);
    out.high = true;
    out.v = value;
    U.ctr(1, 1, 0, nullptr, W.c1.data());
    std::fill(out.d1, out.d1 + mk, 0.0);
    for (int r = 0; r < P; ++r) {
        if (W.c1[r] == 0.0) continue;
        for (int k = 0; k < mk; ++k) out.d1[k] += W.c1[r] * E[r]->d1[k];
    }
    compose2(U, E, W, out.P2);
    for (int d = 0; d < 3; ++d) compose3(U, E, W, d, out.P3[d]);
    for (int pr = 0; pr < 3; ++pr) compose4(U, E, W, pr, out.P4[pr]);
    compose5(U, E, W, out.P5);
    // symmetric by construction, made so exactly
    double* mats[8] = {out.P2, out.P3[0], out.P3[1], out.P3[2], out.P4[0],
                       out.P4[1], out.P4[2], out.P5};
    for (double* h : mats) {
        for (int k = 0; k < mk; ++k) {
            for (int j = k + 1; j < mk; ++j) {
                double s = 0.5 * (h[(std::size_t) k * mk + j] + h[(std::size_t) j * mk + k]);
                h[(std::size_t) k * mk + j] = s;
                h[(std::size_t) j * mk + k] = s;
            }
        }
    }
    out.finish(dirs);
}

// ---- the outer derivatives ---------------------------------------------------

// the model's arrays: u_{r_1..r_k} of the driving score at a row, in
// columns lofs + kmult * (index over the tuples); each array is symmetric
// in its tuple, so the order of the indices does not matter
struct ArrayGet {
    const double* a[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
    int n = 0, P = 0, row = 0;
    std::size_t lofs = 0, kmult = 1;
    double operator()(int k, const int* idx) const {
        std::size_t col = 0;
        for (int i = 0; i < k; ++i) col = col * P + idx[i];
        return a[k - 1][(lofs + kmult * col) * (std::size_t) n + row];
    }
};

// The Levinson-Durbin map's derivatives of orders one to five as full
// tensors over the partial autocorrelations, by its recursion:
// xi^(k)_i = xi^(k-1)_i - rho_k xi^(k-1)_{k-i}, so a derivative of order m
// is the same combination of the derivatives of order m, less the
// derivatives of order m - 1 of xi^(k-1)_{k-i} wherever one of its indices is
// k. T[m][i] has q^m entries, the index (s_1, .., s_m) at sum s_t q^(t-1).
struct Levinson5 {
    int q = 0;
    dvec phi;
    std::vector<dvec> T[6];          // T[m][i], m = 1..5
    std::vector<dvec> nT[6];
    dvec nphi;
    std::vector<int> pw;
    void alloc(int q_) {
        q = q_;
        pw.assign(7, 1);
        for (int m = 1; m <= 6; ++m) pw[m] = pw[m - 1] * q;
        phi.assign(q, 0.0); nphi.assign(q, 0.0);
        for (int m = 1; m <= 5; ++m) {
            T[m].assign(q, dvec(pw[m], 0.0));
            nT[m].assign(q, dvec(pw[m], 0.0));
        }
    }
    void run(const double* rho) {
        int qq = q;
        for (int m = 1; m <= 5; ++m) for (auto& t : T[m]) std::fill(t.begin(), t.end(), 0.0);
        for (int k = 1; k <= qq; ++k) {
            double rk = rho[k - 1];
            for (int m = 1; m <= 5; ++m) for (int i = 0; i < k; ++i) std::fill(nT[m][i].begin(), nT[m][i].end(), 0.0);
            nphi[k - 1] = rk;
            nT[1][k - 1][k - 1] = 1.0;
            for (int i = 0; i < k - 1; ++i) {
                int r = k - 2 - i;
                nphi[i] = phi[i] - rk * phi[r];
                for (int m = 1; m <= 5; ++m) {
                    dvec& out = nT[m][i];
                    const dvec& Ti = T[m][i];
                    const dvec& Tr = T[m][r];
                    for (int c = 0; c < pw[m]; ++c) out[c] = Ti[c] - rk * Tr[c];
                    // the slots that hit rho_k: the order below of xi_r
                    for (int c = 0; c < pw[m]; ++c) {
                        int rem = c;
                        for (int t = 0; t < m; ++t) {
                            int st = rem % qq;
                            rem /= qq;
                            if (st != k - 1) continue;
                            if (m == 1) { out[c] -= phi[r]; continue; }
                            // the index without slot t
                            int lo = c % pw[t];
                            int hi = c / pw[t + 1];
                            int cm = lo + hi * pw[t];
                            out[c] -= T[m - 1][r][cm];
                        }
                    }
                }
            }
            for (int i = 0; i < k; ++i) phi[i] = nphi[i];
            for (int m = 1; m <= 5; ++m) for (int i = 0; i < k; ++i) T[m][i].swap(nT[m][i]);
        }
    }
};

struct LevGet {
    const Levinson5* L;
    int j;
    double operator()(int k, const int* idx) const {
        int c = 0;
        for (int t = k - 1; t >= 0; --t) c = c * L->q + idx[t];
        return L->T[k][j][c];
    }
};

} // namespace mt7f

#endif
