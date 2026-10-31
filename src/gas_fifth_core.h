// The fifth order of the score-driven recursion: the pieces that the three
// compiled kernels (gas_fifth_gen.cpp, gas_fifth_joint.cpp and
// gas_fifth_joint_sub.cpp) share. Everything here is written out by order:
// a value of the recursion carries its derivatives along every subset of the
// three directions (v, w, x), each as an mk x mk matrix with the two free
// slots, and the product of two values is Leibniz's rule over those subsets,
// the sixteen-term rule of the fourth order (gas_curvature_gen.cpp, prod4)
// extended to thirty-two terms. A subset is a bit mask: bit 0 is v, bit 1 is
// w, bit 2 is x.
#ifndef MT7_GAS_FIFTH_CORE_H
#define MT7_GAS_FIFTH_CORE_H

#include <vector>
#include <cstddef>
#include <algorithm>

namespace mt7f {

typedef std::vector<double> dvec;

inline void outer_add(double* h, double c, const double* a, const double* b,
                      int mk) {
    if (c == 0.0) return;
    for (int k = 0; k < mk; ++k) {
        double bk = c * b[k];
        if (bk == 0.0) continue;
        for (int j = 0; j < mk; ++j) h[(std::size_t) k * mk + j] += a[j] * bk;
    }
}
inline void sym_outer_add(double* h, const double* a, const double* b, int mk) {
    for (int k = 0; k < mk; ++k) {
        double ak = a[k], bk = b[k];
        for (int j = 0; j < mk; ++j) {
            h[(std::size_t) k * mk + j] += a[j] * bk + b[j] * ak;
        }
    }
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
inline int top_bit(int S) { return S >= 4 ? 2 : (S >= 2 ? 1 : 0); }
inline int nbits(int S) { return (S & 1) + ((S >> 1) & 1) + ((S >> 2) & 1); }

// A value with its derivatives: v, d1 and the matrix along each subset; the
// rows (one free slot) and the scalars (none) along each subset are derived
// once, by contracting the matrix of the subset less its highest direction
// against that direction.
struct J5 {
    int mk = 0;
    bool linear = false;   // only v and d1: every matrix, row above the first and scalar above the value are zero
    double v = 0.0;
    dvec d1;
    dvec M;          // 8 blocks of mk x mk
    dvec R;          // 8 rows of mk
    double s[8];
    double* mat(int S) { return &M[(std::size_t) S * mk * mk]; }
    const double* mat(int S) const { return &M[(std::size_t) S * mk * mk]; }
    const double* row(int S) const { return &R[(std::size_t) S * mk]; }
    void init(int mk_) {
        mk = mk_;
        d1.assign(mk, 0.0);
        M.assign((std::size_t) 8 * mk * mk, 0.0);
        R.assign((std::size_t) 8 * mk, 0.0);
        std::fill(s, s + 8, 0.0);
        v = 0.0;
    }
    void clear() {
        v = 0.0;
        std::fill(d1.begin(), d1.end(), 0.0);
        std::fill(M.begin(), M.end(), 0.0);
    }
    // the rows and scalars from the matrices and the directions
    void finish(const double* const* dirs) {
        std::copy(d1.begin(), d1.end(), R.begin());
        s[0] = v;
        for (int S = 1; S < 8; ++S) {
            int h = top_bit(S);
            mv(mat(S ^ (1 << h)), dirs[h], &R[(std::size_t) S * mk], mk);
        }
        for (int S = 1; S < 8; ++S) {
            int h = top_bit(S);
            s[S] = dt(row(S ^ (1 << h)), dirs[h], mk);
        }
    }
};

// a predictor of the first order: value, first derivative and the scalars
// along one direction; every other entry is zero
inline void linear_jet(double v0, const double* d1, const double* const* dirs,
                       J5& J) {
    int mk = J.mk;
    J.linear = true;
    J.v = v0;
    std::copy(d1, d1 + mk, J.d1.begin());
    std::fill(J.R.begin(), J.R.end(), 0.0);
    std::copy(d1, d1 + mk, J.R.begin());
    std::fill(J.s, J.s + 8, 0.0);
    J.s[0] = v0;
    for (int d = 0; d < 3; ++d) J.s[1 << d] = dt(d1, dirs[d], mk);
}

// out += the derivative of A B along the subset S, two free slots: the
// Leibniz rule over the subsets T of S
inline void prod_add(const J5& A, const J5& B, int S, double* out) {
    int mk = A.mk;
    std::size_t mk2 = (std::size_t) mk * mk;
    for (int T = S;; T = (T - 1) & S) {
        int U = S ^ T;
        double a0 = A.s[T], b0 = B.s[U];
        const double* a2 = A.mat(T);
        const double* b2 = B.mat(U);
        for (std::size_t e = 0; e < mk2; ++e) out[e] += a2[e] * b0 + a0 * b2[e];
        sym_outer_add(out, A.row(T), B.row(U), mk);
        if (T == 0) break;
    }
}

// C += the full product of A and B: value, first derivative and every
// subset's matrix (the rows and scalars are left to finish())
inline void prod_into(const J5& A, const J5& B, J5& C) {
    C.v += A.v * B.v;
    for (int k = 0; k < A.mk; ++k) C.d1[k] += A.v * B.d1[k] + B.v * A.d1[k];
    for (int S = 0; S < 8; ++S) prod_add(A, B, S, C.mat(S));
}

// the jet of h(z'u), rank one: value k0, first derivative k1 z, and along
// the subset S, k_{2+|S|} prod_{d in S} (z.d) z z'
inline void rank_one(const double* K, const double* z, const double* zd,
                     J5& J) {
    int mk = J.mk;
    J.v = K[0];
    for (int k = 0; k < mk; ++k) J.d1[k] = K[1] * z[k];
    for (int S = 0; S < 8; ++S) {
        double c = K[2 + nbits(S)];
        for (int d = 0; d < 3; ++d) if (S & (1 << d)) c *= zd[d];
        double* h = J.mat(S);
        std::fill(h, h + (std::size_t) mk * mk, 0.0);
        outer_add(h, c, z, z, mk);
    }
}

// the starting level from the fixed point f0 = omega + S f0: along each
// subset, by increasing size, (1 - S)^(-1) times the level's derivative plus
// the product rule of S f0 without its term S * (f0's derivative there)
inline void start_level(const J5& om, const J5& Ssum, J5& f0,
                        const double* const* dirs) {
    int mk = om.mk;
    std::size_t mk2 = (std::size_t) mk * mk;
    double cst = 1.0 / (1.0 - Ssum.v);
    f0.clear();
    f0.v = om.v * cst;
    for (int k = 0; k < mk; ++k) f0.d1[k] = cst * (om.d1[k] + Ssum.d1[k] * f0.v);
    static const int order[8] = {0, 1, 2, 4, 3, 5, 6, 7};
    for (int oi = 0; oi < 8; ++oi) {
        int S = order[oi];
        f0.finish(dirs);          // the lower subsets are final
        dvec acc(mk2, 0.0);
        prod_add(Ssum, f0, S, acc.data());
        const double* o = om.mat(S);
        double* h = f0.mat(S);
        // f0's own matrix at S is still zero, so prod_add left out S * f0_S
        for (std::size_t e = 0; e < mk2; ++e) h[e] = cst * (o[e] + acc[e]);
    }
    f0.finish(dirs);
}

// The model's derivative arrays at one row: u_{r_1 .. r_k} for k = 1..5,
// one column per index tuple (first index most significant), read for the
// driving score `l` of a joint filter or the single score of a term.
struct Arrays {
    int npar = 0;
    int n = 0;
    const double* a[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
};

// the contraction of u^{(k)} at one row against vectors over the
// equations, leaving 0, 1 or 2 indices free (the free ones first)
struct Contract {
    const Arrays* A;
    int row;
    // the column of a tuple is lofs + kmult * (its index over the tuples):
    // lofs = 0, kmult = 1 for a term's arrays (H, D3, ...), and the driving
    // score l with kmult = K for a joint filter's (U1, U2, ...), whose R
    // arrays are n x K x P x ... in column-major order. Every array is
    // symmetric in its tuple, so the order of the indices does not matter.
    std::size_t lofs = 0;
    std::size_t kmult = 1;
    double get(int k, const int* idx) const {
        std::size_t col = 0;
        for (int i = 0; i < k; ++i) col = col * A->npar + idx[i];
        return A->a[k - 1][(lofs + kmult * col) * (std::size_t) A->n + row];
    }
    // sum over the contracted indices; nfree leading indices kept
    void run(int nfree, int nfixed, const double* const* vecs, double* out) const {
        int P = A->npar;
        int k = nfree + nfixed;
        int nout = nfree == 0 ? 1 : (nfree == 1 ? P : P * P);
        for (int o = 0; o < nout; ++o) out[o] = 0.0;
        int idx[5];
        int total = 1;
        for (int i = 0; i < k; ++i) total *= P;
        for (int c = 0; c < total; ++c) {
            int rem = c;
            for (int i = k - 1; i >= 0; --i) { idx[i] = rem % P; rem /= P; }
            double w = 1.0;
            for (int i = 0; i < nfixed; ++i) {
                w *= vecs[i][idx[nfree + i]];
                if (w == 0.0) break;
            }
            if (w == 0.0) continue;
            double val = get(k, idx);
            int o = nfree == 0 ? 0 : (nfree == 1 ? idx[0] : idx[0] + P * idx[1]);
            out[o] += w * val;
        }
    }
};

// the partitions of a bit mask R of at most three directions into blocks
inline void partitions(int R, std::vector<std::vector<int>>& out) {
    out.clear();
    int el[3], ne = 0;
    for (int d = 0; d < 3; ++d) if (R & (1 << d)) el[ne++] = 1 << d;
    if (ne == 0) { out.push_back({}); return; }
    if (ne == 1) { out.push_back({el[0]}); return; }
    if (ne == 2) {
        out.push_back({el[0], el[1]});
        out.push_back({el[0] | el[1]});
        return;
    }
    out.push_back({el[0], el[1], el[2]});
    out.push_back({el[0] | el[1], el[2]});
    out.push_back({el[0] | el[2], el[1]});
    out.push_back({el[1] | el[2], el[0]});
    out.push_back({el[0] | el[1] | el[2]});
}

// The driving score's derivative along the subset S with two free slots,
// written out by Faa di Bruno's formula. E holds the predictors' jets (only
// the filtered equations have derivatives above the first); g[r] is the
// vector of (V_r . d) over the directions d. The free slots either share a
// block (a matrix of one predictor, times the coefficient of the remaining
// directions' partitions), or sit in two blocks (two rows of predictors,
// with a matrix of coefficients over their equations).
inline void score_mat(const Contract& C, const std::vector<const J5*>& E,
                      int S, double* out, int mk) {
    int P = (int) E.size();
    std::size_t mk2 = (std::size_t) mk * mk;
    std::fill(out, out + mk2, 0.0);
    std::vector<std::vector<int>> parts;
    dvec blockvec((std::size_t) 4 * P), coef(P * P), tmp(P * P), acol(mk);
    std::vector<const double*> vecs(4);
    // a block of directions B, no free slot: a vector over the equations
    auto block_vec = [&](int B, double* u) {
        for (int r = 0; r < P; ++r) u[r] = E[r] ? E[r]->s[B] : 0.0;
    };
    // the coefficient over the first nfree indices, summed over the
    // partitions of R
    auto coef_of = [&](int nfree, int R, double* outc) {
        int nout = nfree == 1 ? P : P * P;
        for (int o = 0; o < nout; ++o) outc[o] = 0.0;
        partitions(R, parts);
        for (const auto& pi : parts) {
            int nb = (int) pi.size();
            for (int b = 0; b < nb; ++b) {
                block_vec(pi[b], &blockvec[(std::size_t) b * P]);
                vecs[b] = &blockvec[(std::size_t) b * P];
            }
            C.run(nfree, nb, vecs.data(), tmp.data());
            for (int o = 0; o < nout; ++o) outc[o] += tmp[o];
        }
    };
    // case one: both free slots in one block with the directions B0
    for (int B0 = S;; B0 = (B0 - 1) & S) {
        int R = S ^ B0;
        coef_of(1, R, coef.data());
        for (int r = 0; r < P; ++r) {
            if (!E[r] || E[r]->linear || coef[r] == 0.0) continue;
            // a block of two or more slots on a predictor of the first order
            // only is zero, which its zero matrix gives
            const double* X = E[r]->mat(B0);
            for (std::size_t e = 0; e < mk2; ++e) out[e] += coef[r] * X[e];
        }
        if (B0 == 0) break;
    }
    // case two: the free slot a with Ba, b with Bb, the rest R
    for (int Ba = S;; Ba = (Ba - 1) & S) {
        int rest = S ^ Ba;
        for (int Bb = rest;; Bb = (Bb - 1) & rest) {
            int R = rest ^ Bb;
            coef_of(2, R, coef.data());
            // sum_rb (sum_ra c[ra, rb] row_ra(Ba)) row_rb(Bb)': the rows of a
            // predictor of the first order are zero above Ba = 0
            for (int rb = 0; rb < P; ++rb) {
                if (!E[rb] || (Bb != 0 && E[rb]->linear)) continue;
                std::fill(acol.begin(), acol.end(), 0.0);
                bool any = false;
                for (int ra = 0; ra < P; ++ra) {
                    if (!E[ra] || (Ba != 0 && E[ra]->linear)) continue;
                    double c = coef[ra + P * rb];
                    if (c == 0.0) continue;
                    const double* ro = E[ra]->row(Ba);
                    for (int k = 0; k < mk; ++k) acol[k] += c * ro[k];
                    any = true;
                }
                if (any) outer_add(out, 1.0, acol.data(), E[rb]->row(Bb), mk);
            }
            if (Bb == 0) break;
        }
        if (Ba == 0) break;
    }
    // symmetric by construction; made so exactly
    for (int k = 0; k < mk; ++k) {
        for (int j = k + 1; j < mk; ++j) {
            double a = 0.5 * (out[(std::size_t) k * mk + j] + out[(std::size_t) j * mk + k]);
            out[(std::size_t) k * mk + j] = a;
            out[(std::size_t) j * mk + k] = a;
        }
    }
}

// the score's jet: value given, first derivative sum_r u_r V_r, and every
// subset's matrix
inline void score_jet(const Contract& C, const std::vector<const J5*>& E,
                      double value, J5& out, const double* const* dirs) {
    int P = (int) E.size();
    int mk = out.mk;
    out.v = value;
    dvec u1(P);
    C.run(1, 0, nullptr, u1.data());
    std::fill(out.d1.begin(), out.d1.end(), 0.0);
    for (int r = 0; r < P; ++r) {
        if (!E[r] || u1[r] == 0.0) continue;
        for (int k = 0; k < mk; ++k) out.d1[k] += u1[r] * E[r]->d1[k];
    }
    for (int S = 0; S < 8; ++S) score_mat(C, E, S, out.mat(S), mk);
    out.linear = false;
    out.finish(dirs);
}

} // namespace mt7f

#endif
