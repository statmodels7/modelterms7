#include <Rcpp.h>
#include <vector>
#include "gas_fifth_core.h"
using namespace Rcpp;

// The fifth order of a score-driven term's recursion, compiled, for the
// route without subformulas and the route with them. The level and the
// loadings are read per row as h(z'u), with z the row of their design among
// the group's active unknowns (a unit vector for a parameter without a
// subformula) and K the link's derivatives of orders zero to five there.
// The autoregressive coefficients are constant over the rows and come with
// their derivatives along every subset of the directions, from R. The
// driving score's derivatives are read from the arrays H, D3, D4, D5 and D6.
//
// The state of a row is the level with its derivatives along every subset of
// the three directions, the score likewise, in ring buffers of max(p, q) + 1
// slots.

// [[Rcpp::export]]
List gas_fifth_gen_cpp(NumericVector eta, List groups, int p, int q,
                       NumericMatrix dirs_full, int ap, NumericVector s_at,
                       Nullable<Function> score_fn, NumericVector g,
                       List arrays, List Vs, NumericMatrix seed) {
    using namespace mt7f;
    int n = eta.size();
    int m_full = seed.ncol();
    int P = Vs.size();
    NumericMatrix D(n, m_full), dom(n, m_full);
    std::vector<NumericMatrix> dphi, dpsi;
    for (int d = 0; d < 3; ++d) {
        dphi.push_back(NumericMatrix(n, m_full));
        dpsi.push_back(NumericMatrix(n, m_full));
    }
    NumericMatrix W(m_full, m_full);
    Arrays Ar;
    Ar.npar = P;
    Ar.n = n;
    std::vector<NumericMatrix> keep_arr;
    for (int k = 0; k < 5; ++k) {
        NumericMatrix x = arrays[k];
        keep_arr.push_back(x);
        Ar.a[k] = keep_arr.back().begin();
    }
    std::vector<NumericMatrix> Vs_m;
    for (int r = 0; r < P; ++r) Vs_m.push_back(as<NumericMatrix>(Vs[r]));
    bool use_fn = s_at.size() == 0;
    Function* fn = nullptr;
    Function fn_obj("identity");
    if (use_fn) {
        fn_obj = Function(score_fn.get());
        fn = &fn_obj;
    }
    const int pairs[3] = {3, 5, 6};

    for (int gi = 0; gi < groups.size(); ++gi) {
        List G = groups[gi];
        IntegerVector rows = G["rows"];
        IntegerVector act = G["act"];
        int m = rows.size(), mk = act.size();
        std::size_t mk2 = (std::size_t) mk * mk;
        // the directions on the active set
        dvec dv((std::size_t) 3 * mk);
        for (int d = 0; d < 3; ++d) {
            for (int k = 0; k < mk; ++k) dv[(std::size_t) d * mk + k] = dirs_full(act[k] - 1, d);
        }
        const double* dirs[3] = {&dv[0], &dv[mk], &dv[2 * mk]};
        // the values read per row: Z is m x mk, K is m x 6
        List om = G["omega"];
        List al = G["a"];
        NumericMatrix omZ = om["Z"], omK = om["K"];
        std::vector<NumericMatrix> aZ, aK;
        for (int i = 0; i < p; ++i) {
            List ai = al[i];
            aZ.push_back(as<NumericMatrix>(ai["Z"]));
            aK.push_back(as<NumericMatrix>(ai["K"]));
        }
        // the constant autoregressive coefficients' jets
        std::vector<J5> bj(q);
        J5 Ssum;
        Ssum.init(mk);
        if (q > 0) {
            NumericVector bval = G["b"];
            NumericMatrix bd1 = G["b_d1"];       // mk x q
            NumericVector bmats = G["b_mats"];   // mk x mk x 8 x q
            for (int j = 0; j < q; ++j) {
                bj[j].init(mk);
                bj[j].v = bval[j];
                for (int k = 0; k < mk; ++k) bj[j].d1[k] = bd1(k, j);
                std::copy(bmats.begin() + (std::size_t) j * 8 * mk2,
                          bmats.begin() + (std::size_t)(j + 1) * 8 * mk2,
                          bj[j].M.begin());
                bj[j].finish(dirs);
                Ssum.v += bj[j].v;
                for (int k = 0; k < mk; ++k) Ssum.d1[k] += bj[j].d1[k];
                for (std::size_t e = 0; e < 8 * mk2; ++e) Ssum.M[e] += bj[j].M[e];
            }
            Ssum.finish(dirs);
        }
        bool direct = as<bool>(G["direct"]);
        double d0 = as<double>(G["d0"]);

        auto value_jet = [&](const NumericMatrix& Z, const NumericMatrix& K,
                             int t, J5& J) {
            dvec z(mk), Kt(6), zd(3);
            for (int k = 0; k < mk; ++k) z[k] = Z(t, k);
            for (int o = 0; o < 6; ++o) Kt[o] = K(t, o);
            for (int d = 0; d < 3; ++d) zd[d] = dt(z.data(), dirs[d], mk);
            rank_one(Kt.data(), z.data(), zd.data(), J);
            J.finish(dirs);
        };

        int L = std::max(p, q) + 1;
        std::vector<J5> fj(L), sj(L);
        for (int i = 0; i < L; ++i) { fj[i].init(mk); sj[i].init(mk); }
        J5 f0, omj, aj, ft, ej, zero_j;
        f0.init(mk); omj.init(mk); aj.init(mk); ft.init(mk); ej.init(mk);
        zero_j.init(mk);
        zero_j.finish(dirs);
        if (direct) {
            f0.v = d0;
            for (int k = 0; k < mk; ++k) f0.d1[k] = -seed(rows[0] - 1, act[k] - 1);
            f0.finish(dirs);
        } else {
            value_jet(omZ, omK, 0, omj);
            start_level(omj, Ssum, f0, dirs);
        }
        std::vector<J5> lin(P);
        for (int r = 0; r < P; ++r) lin[r].init(mk);
        std::vector<const J5*> E(P);
        dvec Wl(mk2, 0.0), vrow(mk);

        for (int t = 0; t < m; ++t) {
            int row = rows[t] - 1;
            std::size_t ts = (std::size_t)(t % L);
            if (direct && t == 0) {
                ft = f0;
            } else {
                value_jet(omZ, omK, t, omj);
                ft = omj;
                for (int i = 1; i <= p; ++i) {
                    int lag = t - i;
                    value_jet(aZ[i - 1], aK[i - 1], t, aj);
                    prod_into(aj, lag >= 0 ? sj[lag % L] : zero_j, ft);
                }
                for (int j = 1; j <= q; ++j) {
                    int lag = t - j;
                    prod_into(bj[j - 1], lag >= 0 ? fj[lag % L] : f0, ft);
                }
                ft.finish(dirs);
            }
            fj[ts] = ft;
            ej = ft;
            ej.v += eta[row];
            for (int k = 0; k < mk; ++k) ej.d1[k] += seed(row, act[k] - 1);
            ej.finish(dirs);
            for (int k = 0; k < mk; ++k) D(row, act[k] - 1) = ej.d1[k];
            for (int d = 0; d < 3; ++d) {
                const double* r1 = ej.row(1 << d);
                const double* r2 = ej.row(pairs[d]);
                for (int k = 0; k < mk; ++k) {
                    dphi[d](row, act[k] - 1) = r1[k];
                    dpsi[d](row, act[k] - 1) = r2[k];
                }
            }
            const double* r3 = ej.row(7);
            for (int k = 0; k < mk; ++k) dom(row, act[k] - 1) = r3[k];
            const double* top = ft.mat(7);
            for (std::size_t e = 0; e < mk2; ++e) Wl[e] += g[row] * top[e];

            // the driving score's jet at this row
            for (int r = 0; r < P; ++r) {
                if (r == ap - 1) { E[r] = &ej; continue; }
                for (int k = 0; k < mk; ++k) vrow[k] = Vs_m[r](row, act[k] - 1);
                linear_jet(0.0, vrow.data(), dirs, lin[r]);
                E[r] = &lin[r];
            }
            double sv = use_fn ? as<double>((*fn)(ej.v, row + 1)) : s_at[row];
            Contract C;
            C.A = &Ar;
            C.row = row;
            score_jet(C, E, sv, sj[ts], dirs);
        }
        for (int j = 0; j < mk; ++j) {
            for (int k = 0; k < mk; ++k) {
                W(act[j] - 1, act[k] - 1) += Wl[(std::size_t) k * mk + j];
            }
        }
    }
    // symmetric to a rounding, made so exactly
    for (int j = 0; j < m_full; ++j) {
        for (int k = j + 1; k < m_full; ++k) {
            double a = 0.5 * (W(j, k) + W(k, j));
            W(j, k) = a;
            W(k, j) = a;
        }
    }
    List dphl(3), dpsl(3);
    for (int d = 0; d < 3; ++d) { dphl[d] = dphi[d]; dpsl[d] = dpsi[d]; }
    return List::create(_["jacobian"] = D, _["dphi"] = dphl, _["dpsi"] = dpsl,
                        _["domega"] = dom, _["curvature"] = W);
}
