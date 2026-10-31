#include <Rcpp.h>
#include <vector>
#include "gas_fifth_core.h"
using namespace Rcpp;

// The fifth order of a score-driven term's recursion, compiled, for the
// route without subformulas and the route with them. Every value of the
// recursion (the level, the loadings, the partial autocorrelations) is read
// per row as h(z'u), with z its design row among the group's active unknowns
// (a unit vector for a parameter without a subformula) and K its link's
// derivatives of orders zero to five there. The autoregressive coefficients
// are the Levinson-Durbin map of the partial autocorrelations, composed by
// order (compose2 .. compose5 of gas_fifth_core.h) once per group where the
// partial autocorrelations are constant and once per row where one carries a
// subformula. The driving score is composed the same way from the arrays H,
// D3, D4, D5 and D6. All storage is allocated once per group.

namespace {

struct RowVal {
    NumericMatrix Z, K;
};

inline void value_jet(const RowVal& V, int t, int mk, mt7f::dvec& z,
                      mt7f::dvec& zd, double* Kt, const double* const* dirs,
                      mt7f::Jet5& J) {
    for (int k = 0; k < mk; ++k) z[k] = V.Z(t, k);
    for (int o = 0; o < 6; ++o) Kt[o] = V.K(t, o);
    for (int d = 0; d < 3; ++d) zd[d] = mt7f::dt(z.data(), dirs[d], mk);
    mt7f::rank_one(Kt, z.data(), zd.data(), J);
    J.finish(dirs);
}

} // namespace

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
    std::vector<NumericMatrix> keep;
    ArrayGet AG;
    AG.n = n;
    AG.P = P;
    for (int k = 0; k < 5; ++k) {
        keep.push_back(as<NumericMatrix>(arrays[k]));
        AG.a[k] = keep.back().begin();
    }
    std::vector<NumericMatrix> Vs_m;
    for (int r = 0; r < P; ++r) Vs_m.push_back(as<NumericMatrix>(Vs[r]));
    bool use_fn = s_at.size() == 0;
    Function fn_obj("identity");
    if (use_fn) fn_obj = Function(score_fn.get());
    Outer<ArrayGet> US;
    US.get = &AG;
    US.P = P;

    for (int gi = 0; gi < groups.size(); ++gi) {
        List G = groups[gi];
        IntegerVector rows = G["rows"];
        IntegerVector act = G["act"];
        int m = rows.size(), mk = act.size();
        std::size_t mk2 = (std::size_t) mk * mk;
        dvec dv((std::size_t) 3 * mk);
        for (int d = 0; d < 3; ++d) {
            for (int k = 0; k < mk; ++k) dv[(std::size_t) d * mk + k] = dirs_full(act[k] - 1, d);
        }
        const double* dirs[3] = {&dv[0], &dv[mk], &dv[2 * mk]};
        RowVal om;
        { List x = G["omega"]; om.Z = as<NumericMatrix>(x["Z"]); om.K = as<NumericMatrix>(x["K"]); }
        std::vector<RowVal> al(p), pa(q);
        { List x = G["a"]; for (int i = 0; i < p; ++i) { List y = x[i]; al[i].Z = as<NumericMatrix>(y["Z"]); al[i].K = as<NumericMatrix>(y["K"]); } }
        { List x = G["pa"]; for (int j = 0; j < q; ++j) { List y = x[j]; pa[j].Z = as<NumericMatrix>(y["Z"]); pa[j].K = as<NumericMatrix>(y["K"]); } }
        bool varying_b = as<bool>(G["varying_b"]);
        bool direct = as<bool>(G["direct"]);
        double d0 = as<double>(G["d0"]);

        // the storage of the group, allocated once
        dvec z(mk), zd(3), Kt(6), rho(q), vrow(mk), Wl(mk2, 0.0);
        int L = std::max(p, q) + 1;
        std::vector<Jet5> fj(L), sj(L), rj(q), bj(q), lin(P);
        for (auto& x : fj) x.alloc(mk);
        for (auto& x : sj) x.alloc(mk);
        for (auto& x : rj) x.alloc(mk);
        for (auto& x : bj) x.alloc(mk);
        for (auto& x : lin) x.alloc(mk);
        Jet5 omj, aj, ft, ej, f0, Ssum, zero_j;
        omj.alloc(mk); aj.alloc(mk); ft.alloc(mk); ej.alloc(mk); f0.alloc(mk);
        Ssum.alloc(mk); zero_j.alloc(mk);
        zero_j.finish(dirs);
        CompWork CW, CB;
        CW.alloc(P, mk);
        CB.alloc(std::max(q, 1), mk);
        Levinson5 Lv;
        Lv.alloc(std::max(q, 1));
        std::vector<const Jet5*> Eb(q), E(P);
        for (int j = 0; j < q; ++j) Eb[j] = &rj[j];

        // the autoregressive coefficients at row t, composed from the
        // partial autocorrelations
        auto b_at = [&](int t) {
            for (int k = 0; k < q; ++k) {
                value_jet(pa[k], t, mk, z, zd, Kt.data(), dirs, rj[k]);
                rho[k] = rj[k].v;
            }
            Lv.run(rho.data());
            Ssum.clear();
            for (int j = 0; j < q; ++j) {
                LevGet LG;
                LG.L = &Lv;
                LG.j = j;
                Outer<LevGet> UB;
                UB.get = &LG;
                UB.P = q;
                compose(UB, Eb, Lv.phi[j], CB, bj[j], dirs);
                Ssum.v += bj[j].v;
                for (int k = 0; k < mk; ++k) Ssum.d1[k] += bj[j].d1[k];
                for (std::size_t e = 0; e < 8 * mk2; ++e) Ssum.buf[e] += bj[j].buf[e];
            }
            Ssum.finish(dirs);
        };
        if (q > 0) b_at(0);
        if (direct) {
            f0.clear();
            f0.v = d0;
            for (int k = 0; k < mk; ++k) f0.d1[k] = -seed(rows[0] - 1, act[k] - 1);
            f0.finish(dirs);
        } else {
            value_jet(om, 0, mk, z, zd, Kt.data(), dirs, omj);
            start_level(omj, Ssum, f0, dirs);
        }

        for (int t = 0; t < m; ++t) {
            int row = rows[t] - 1;
            std::size_t ts = (std::size_t)(t % L);
            if (direct && t == 0) {
                ft = f0;
            } else {
                if (q > 0 && varying_b && t > 0) b_at(t);
                value_jet(om, t, mk, z, zd, Kt.data(), dirs, ft);
                for (int i = 1; i <= p; ++i) {
                    int lag = t - i;
                    value_jet(al[i - 1], t, mk, z, zd, Kt.data(), dirs, aj);
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
                const double* r1 = ej.r2[d];
                const double* r2 = ej.r3[d];
                for (int k = 0; k < mk; ++k) {
                    dphi[d](row, act[k] - 1) = r1[k];
                    dpsi[d](row, act[k] - 1) = r2[k];
                }
            }
            for (int k = 0; k < mk; ++k) dom(row, act[k] - 1) = ej.r4[k];
            for (std::size_t e = 0; e < mk2; ++e) Wl[e] += g[row] * ft.P5[e];

            // the driving score's jet at this row
            for (int r = 0; r < P; ++r) {
                if (r == ap - 1) { E[r] = &ej; continue; }
                Jet5& Lr = lin[r];
                Lr.high = false;
                Lr.v = 0.0;
                for (int k = 0; k < mk; ++k) Lr.d1[k] = Vs_m[r](row, act[k] - 1);
                Lr.finish(dirs);
                E[r] = &Lr;
            }
            double sv = use_fn ? as<double>(fn_obj(ej.v, row + 1)) : s_at[row];
            AG.row = row;
            compose(US, E, sv, CW, sj[ts], dirs);
        }
        for (int j = 0; j < mk; ++j) {
            for (int k = 0; k < mk; ++k) {
                W(act[j] - 1, act[k] - 1) += Wl[(std::size_t) k * mk + j];
            }
        }
    }
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
