#include <Rcpp.h>
#include <vector>
#include "gas_fifth_core.h"
using namespace Rcpp;

// The fifth order of a joint score-driven filter, compiled, for members with
// and without subformulas. Each member's level, loadings, partial
// autocorrelations and cross loadings are read per row as h(z'u) through
// their design rows on the group's active unknowns and their links'
// derivatives of orders zero to five; the autoregressive coefficients are
// the Levinson-Durbin map of the partial autocorrelations, composed by order
// once per group, or once per row where a partial autocorrelation carries a
// subformula. A driver (a parameter joined without a recursion of its own)
// has a level of zero. The driving scores are composed by order from U1 ..
// U5, n x K x P^k arrays. All storage is allocated once per group.

namespace {

struct RowVal {
    NumericMatrix Z, K;
};

struct Member {
    bool driver = false, direct = false, varying_b = false;
    int p = 0, q = 0, eq = 0;
    double d0 = 0.0;
    RowVal om;
    std::vector<RowVal> al, pa, xs;
    std::vector<int> xpar, xlag;
    std::vector<mt7f::Jet5> rj, bj;
    mt7f::Jet5 Ssum;
    mt7f::Levinson5 Lv;
    mt7f::CompWork CB;
    std::vector<const mt7f::Jet5*> Eb;
};

RowVal read_rv(List x) {
    RowVal V;
    V.Z = as<NumericMatrix>(x["Z"]);
    V.K = as<NumericMatrix>(x["K"]);
    return V;
}

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
List gas_fifth_joint_cpp(List groups, int K, NumericMatrix dirs_full,
                         List seeds, IntegerVector fp, NumericMatrix g,
                         NumericMatrix u_at, List arrays) {
    using namespace mt7f;
    int P = seeds.size();
    std::vector<NumericMatrix> seed;
    for (int a = 0; a < P; ++a) seed.push_back(as<NumericMatrix>(seeds[a]));
    int n = seed[0].nrow();
    int m_full = seed[0].ncol();
    std::vector<NumericVector> keep;
    ArrayGet AG;
    AG.n = n;
    AG.P = P;
    AG.kmult = (std::size_t) K;
    for (int k = 0; k < 5; ++k) {
        keep.push_back(as<NumericVector>(arrays[k]));
        AG.a[k] = keep.back().begin();
    }
    Outer<ArrayGet> US;
    US.get = &AG;
    US.P = P;
    std::vector<NumericMatrix> J_(K), Dm(K);
    std::vector<std::vector<NumericMatrix>> Ph(K), Ps(K);
    for (int k = 0; k < K; ++k) {
        J_[k] = NumericMatrix(n, m_full);
        Dm[k] = NumericMatrix(n, m_full);
        for (int d = 0; d < 3; ++d) {
            Ph[k].push_back(NumericMatrix(n, m_full));
            Ps[k].push_back(NumericMatrix(n, m_full));
        }
    }
    NumericMatrix W(m_full, m_full);

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
        dvec z(mk), zd(3), Kt(6), rho, vrow(mk), Wl(mk2, 0.0);
        List mem = G["members"];
        std::vector<Member> M(K);
        int L = 1;
        for (int k = 0; k < K; ++k) {
            List mk_ = mem[k];
            Member& Mb = M[k];
            Mb.eq = fp[k] - 1;
            Mb.driver = as<bool>(mk_["driver"]);
            Mb.direct = as<bool>(mk_["direct"]);
            Mb.d0 = as<double>(mk_["d0"]);
            Mb.Ssum.alloc(mk);
            Mb.Ssum.finish(dirs);
            if (Mb.driver) continue;
            Mb.p = as<int>(mk_["p"]);
            Mb.q = as<int>(mk_["q"]);
            Mb.varying_b = as<bool>(mk_["varying_b"]);
            Mb.om = read_rv(as<List>(mk_["omega"]));
            List al = mk_["a"], pa = mk_["pa"], xl = mk_["x"];
            for (int i = 0; i < Mb.p; ++i) Mb.al.push_back(read_rv(as<List>(al[i])));
            for (int j = 0; j < Mb.q; ++j) Mb.pa.push_back(read_rv(as<List>(pa[j])));
            IntegerVector xpar = mk_["xpar"], xlag = mk_["xlag"];
            for (int h = 0; h < xl.size(); ++h) {
                Mb.xs.push_back(read_rv(as<List>(xl[h])));
                Mb.xpar.push_back(xpar[h] - 1);
                Mb.xlag.push_back(xlag[h]);
                L = std::max(L, xlag[h] + 1);
            }
            L = std::max(L, std::max(Mb.p, Mb.q) + 1);
            Mb.rj.resize(Mb.q);
            Mb.bj.resize(Mb.q);
            for (int j = 0; j < Mb.q; ++j) { Mb.rj[j].alloc(mk); Mb.bj[j].alloc(mk); }
            Mb.Lv.alloc(std::max(Mb.q, 1));
            Mb.CB.alloc(std::max(Mb.q, 1), mk);
            Mb.Eb.resize(Mb.q);
            for (int j = 0; j < Mb.q; ++j) Mb.Eb[j] = &Mb.rj[j];
        }
        int qmax = 1;
        for (int k = 0; k < K; ++k) qmax = std::max(qmax, M[k].q);
        rho.assign(qmax, 0.0);
        auto b_at = [&](Member& Mb, int t) {
            for (int k = 0; k < Mb.q; ++k) {
                value_jet(Mb.pa[k], t, mk, z, zd, Kt.data(), dirs, Mb.rj[k]);
                rho[k] = Mb.rj[k].v;
            }
            Mb.Lv.run(rho.data());
            Mb.Ssum.clear();
            for (int j = 0; j < Mb.q; ++j) {
                LevGet LG;
                LG.L = &Mb.Lv;
                LG.j = j;
                Outer<LevGet> UB;
                UB.get = &LG;
                UB.P = Mb.q;
                compose(UB, Mb.Eb, Mb.Lv.phi[j], Mb.CB, Mb.bj[j], dirs);
                Mb.Ssum.v += Mb.bj[j].v;
                for (int c = 0; c < mk; ++c) Mb.Ssum.d1[c] += Mb.bj[j].d1[c];
                for (std::size_t e = 0; e < 8 * mk2; ++e) Mb.Ssum.buf[e] += Mb.bj[j].buf[e];
            }
            Mb.Ssum.finish(dirs);
        };
        // the storage of the group, allocated once
        Jet5 tmp, zero_j, ft;
        tmp.alloc(mk); zero_j.alloc(mk); ft.alloc(mk);
        zero_j.finish(dirs);
        std::vector<Jet5> f0(K);
        std::vector<std::vector<Jet5>> fj(K, std::vector<Jet5>(L)), uj(K, std::vector<Jet5>(L));
        for (int k = 0; k < K; ++k) {
            f0[k].alloc(mk);
            for (int i = 0; i < L; ++i) { fj[k][i].alloc(mk); uj[k][i].alloc(mk); }
        }
        std::vector<Jet5> Ej(P);
        for (int a = 0; a < P; ++a) Ej[a].alloc(mk);
        std::vector<const Jet5*> E(P);
        CompWork CW;
        CW.alloc(P, mk);
        // the starting levels
        for (int k = 0; k < K; ++k) {
            Member& Mb = M[k];
            if (!Mb.driver && Mb.q > 0) b_at(Mb, 0);
            if (Mb.direct) {
                f0[k].clear();
                f0[k].v = Mb.d0;
                for (int c = 0; c < mk; ++c) f0[k].d1[c] = -seed[Mb.eq](rows[0] - 1, act[c] - 1);
                f0[k].finish(dirs);
            } else if (Mb.driver) {
                f0[k].clear();
                f0[k].finish(dirs);
            } else {
                value_jet(Mb.om, 0, mk, z, zd, Kt.data(), dirs, tmp);
                start_level(tmp, Mb.Ssum, f0[k], dirs);
            }
        }

        for (int t = 0; t < m; ++t) {
            int row = rows[t] - 1;
            std::size_t ts = (std::size_t)(t % L);
            for (int k = 0; k < K; ++k) {
                Member& Mb = M[k];
                if (Mb.direct && t == 0) { fj[k][ts] = f0[k]; continue; }
                if (Mb.driver) { fj[k][ts] = zero_j; continue; }
                if (Mb.q > 0 && Mb.varying_b && t > 0) b_at(Mb, t);
                value_jet(Mb.om, t, mk, z, zd, Kt.data(), dirs, ft);
                for (int i = 1; i <= Mb.p; ++i) {
                    if (t - i < 0) continue;
                    value_jet(Mb.al[i - 1], t, mk, z, zd, Kt.data(), dirs, tmp);
                    prod_into(tmp, uj[k][(t - i) % L], ft);
                }
                for (std::size_t h = 0; h < Mb.xpar.size(); ++h) {
                    if (t - Mb.xlag[h] < 0) continue;
                    value_jet(Mb.xs[h], t, mk, z, zd, Kt.data(), dirs, tmp);
                    prod_into(tmp, uj[Mb.xpar[h]][(t - Mb.xlag[h]) % L], ft);
                }
                for (int j = 1; j <= Mb.q; ++j) {
                    prod_into(Mb.bj[j - 1], t - j >= 0 ? fj[k][(t - j) % L] : f0[k], ft);
                }
                ft.finish(dirs);
                fj[k][ts] = ft;
            }
            // the predictors of every equation at this row
            for (int a = 0; a < P; ++a) {
                Jet5& Ea = Ej[a];
                Ea.high = false;
                Ea.v = 0.0;
                for (int c = 0; c < mk; ++c) Ea.d1[c] = seed[a](row, act[c] - 1);
            }
            for (int k = 0; k < K; ++k) {
                Jet5& Ea = Ej[M[k].eq];
                const Jet5& F = fj[k][ts];
                for (int c = 0; c < mk; ++c) Ea.d1[c] += F.d1[c];
                std::copy(F.buf.begin(), F.buf.begin() + 8 * mk2, Ea.buf.begin());
                Ea.high = true;
            }
            for (int a = 0; a < P; ++a) { Ej[a].finish(dirs); E[a] = &Ej[a]; }
            for (int k = 0; k < K; ++k) {
                const Jet5& Ea = Ej[M[k].eq];
                for (int c = 0; c < mk; ++c) J_[k](row, act[c] - 1) = Ea.d1[c];
                for (int d = 0; d < 3; ++d) {
                    for (int c = 0; c < mk; ++c) {
                        Ph[k][d](row, act[c] - 1) = Ea.r2[d][c];
                        Ps[k][d](row, act[c] - 1) = Ea.r3[d][c];
                    }
                }
                for (int c = 0; c < mk; ++c) Dm[k](row, act[c] - 1) = Ea.r4[c];
                double gk = g(row, k);
                const double* top = fj[k][ts].P5;
                for (std::size_t e = 0; e < mk2; ++e) Wl[e] += gk * top[e];
            }
            AG.row = row;
            for (int l = 0; l < K; ++l) {
                AG.lofs = (std::size_t) l;
                compose(US, E, u_at(row, l), CW, uj[l][ts], dirs);
            }
        }
        for (int j = 0; j < mk; ++j) {
            for (int c = 0; c < mk; ++c) {
                W(act[j] - 1, act[c] - 1) += Wl[(std::size_t) c * mk + j];
            }
        }
    }
    for (int j = 0; j < m_full; ++j) {
        for (int c = j + 1; c < m_full; ++c) {
            double a = 0.5 * (W(j, c) + W(c, j));
            W(j, c) = a;
            W(c, j) = a;
        }
    }
    List jac(K), dom(K), dphi(K), dpsi(K);
    for (int k = 0; k < K; ++k) {
        jac[k] = J_[k];
        dom[k] = Dm[k];
        List a(3), b(3);
        for (int d = 0; d < 3; ++d) { a[d] = Ph[k][d]; b[d] = Ps[k][d]; }
        dphi[k] = a;
        dpsi[k] = b;
    }
    return List::create(_["jacobian"] = jac, _["dphi"] = dphi, _["dpsi"] = dpsi,
                        _["domega"] = dom, _["curvature"] = W);
}
