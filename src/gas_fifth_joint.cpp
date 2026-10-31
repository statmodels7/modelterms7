#include <Rcpp.h>
#include <vector>
#include "gas_fifth_core.h"
using namespace Rcpp;

// The fifth order of a joint score-driven filter, compiled, for members with
// and without subformulas. Each member's level, loadings and cross loadings
// are read per row as h(z'u) through their design rows on the group's active
// unknowns and their links' derivatives of orders zero to five; its
// autoregressive coefficients are constant over the rows and come from R
// with their derivatives along every subset of the directions. A driver
// (a parameter joined without a recursion of its own) has a level of zero.
// The driving scores' derivatives in the predictors of every equation are
// read from U1 .. U5, n x K x P^k arrays.

namespace {

struct Member {
    bool driver = false;
    int p = 0, q = 0, eq = 0;
    NumericMatrix omZ, omK;
    std::vector<NumericMatrix> aZ, aK, xZ, xK;
    std::vector<int> xpar, xlag;
    std::vector<mt7f::J5> bj;
    mt7f::J5 Ssum;
    bool direct = false;
    double d0 = 0.0;
};

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
    Arrays Ar;
    Ar.npar = P;
    Ar.n = n;
    std::vector<NumericVector> keep_arr;
    for (int k = 0; k < 5; ++k) {
        keep_arr.push_back(as<NumericVector>(arrays[k]));
        Ar.a[k] = keep_arr.back().begin();
    }
    List jac(K), dom(K), dphi(K), dpsi(K);
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
    const int pairs[3] = {3, 5, 6};

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
        List mem = G["members"];
        std::vector<Member> M(K);
        for (int k = 0; k < K; ++k) {
            List mk_ = mem[k];
            Member& Mb = M[k];
            Mb.eq = fp[k] - 1;
            Mb.driver = as<bool>(mk_["driver"]);
            Mb.Ssum.init(mk);
            Mb.Ssum.finish(dirs);
            Mb.direct = as<bool>(mk_["direct"]);
            Mb.d0 = as<double>(mk_["d0"]);
            if (Mb.driver) continue;
            Mb.p = as<int>(mk_["p"]);
            Mb.q = as<int>(mk_["q"]);
            List om = mk_["omega"];
            Mb.omZ = as<NumericMatrix>(om["Z"]);
            Mb.omK = as<NumericMatrix>(om["K"]);
            List al = mk_["a"];
            for (int i = 0; i < Mb.p; ++i) {
                List ai = al[i];
                Mb.aZ.push_back(as<NumericMatrix>(ai["Z"]));
                Mb.aK.push_back(as<NumericMatrix>(ai["K"]));
            }
            List xl = mk_["x"];
            IntegerVector xpar = mk_["xpar"], xlag = mk_["xlag"];
            for (int h = 0; h < xl.size(); ++h) {
                List xh = xl[h];
                Mb.xZ.push_back(as<NumericMatrix>(xh["Z"]));
                Mb.xK.push_back(as<NumericMatrix>(xh["K"]));
                Mb.xpar.push_back(xpar[h] - 1);
                Mb.xlag.push_back(xlag[h]);
            }
            if (Mb.q > 0) {
                NumericVector bval = mk_["b"];
                NumericMatrix bd1 = mk_["b_d1"];
                NumericVector bmats = mk_["b_mats"];
                Mb.bj.resize(Mb.q);
                for (int j = 0; j < Mb.q; ++j) {
                    J5& B = Mb.bj[j];
                    B.init(mk);
                    B.v = bval[j];
                    for (int c = 0; c < mk; ++c) B.d1[c] = bd1(c, j);
                    std::copy(bmats.begin() + (std::size_t) j * 8 * mk2,
                              bmats.begin() + (std::size_t)(j + 1) * 8 * mk2,
                              B.M.begin());
                    B.finish(dirs);
                    Mb.Ssum.v += B.v;
                    for (int c = 0; c < mk; ++c) Mb.Ssum.d1[c] += B.d1[c];
                    for (std::size_t e = 0; e < 8 * mk2; ++e) Mb.Ssum.M[e] += B.M[e];
                }
                Mb.Ssum.finish(dirs);
            }
        }
        auto value_jet = [&](const NumericMatrix& Z, const NumericMatrix& Kd,
                             int t, J5& J) {
            dvec z(mk), Kt(6), zd(3);
            for (int c = 0; c < mk; ++c) z[c] = Z(t, c);
            for (int o = 0; o < 6; ++o) Kt[o] = Kd(t, o);
            for (int d = 0; d < 3; ++d) zd[d] = dt(z.data(), dirs[d], mk);
            rank_one(Kt.data(), z.data(), zd.data(), J);
            J.finish(dirs);
        };
        // the starting levels
        std::vector<J5> f0(K);
        J5 tmp, zero_j, ft;
        tmp.init(mk); zero_j.init(mk); ft.init(mk);
        zero_j.finish(dirs);
        for (int k = 0; k < K; ++k) {
            f0[k].init(mk);
            Member& Mb = M[k];
            if (Mb.direct) {
                f0[k].v = Mb.d0;
                for (int c = 0; c < mk; ++c) f0[k].d1[c] = -seed[Mb.eq](rows[0] - 1, act[c] - 1);
                f0[k].finish(dirs);
            } else if (Mb.driver) {
                f0[k].finish(dirs);
            } else {
                value_jet(Mb.omZ, Mb.omK, 0, tmp);
                start_level(tmp, Mb.Ssum, f0[k], dirs);
            }
        }
        int L = 1;
        for (int k = 0; k < K; ++k) {
            L = std::max(L, std::max(M[k].p, M[k].q) + 1);
            for (int lag : M[k].xlag) L = std::max(L, lag + 1);
        }
        std::vector<std::vector<J5>> fj(K, std::vector<J5>(L)), uj(K, std::vector<J5>(L));
        for (int k = 0; k < K; ++k) {
            for (int i = 0; i < L; ++i) { fj[k][i].init(mk); uj[k][i].init(mk); }
        }
        std::vector<J5> Ej(P);
        for (int a = 0; a < P; ++a) Ej[a].init(mk);
        std::vector<const J5*> E(P);
        dvec Wl(mk2, 0.0), vrow(mk);

        for (int t = 0; t < m; ++t) {
            int row = rows[t] - 1;
            std::size_t ts = (std::size_t)(t % L);
            for (int k = 0; k < K; ++k) {
                Member& Mb = M[k];
                if (Mb.direct && t == 0) { fj[k][ts] = f0[k]; continue; }
                if (Mb.driver) { fj[k][ts] = zero_j; continue; }
                value_jet(Mb.omZ, Mb.omK, t, ft);
                for (int i = 1; i <= Mb.p; ++i) {
                    if (t - i < 0) continue;
                    value_jet(Mb.aZ[i - 1], Mb.aK[i - 1], t, tmp);
                    prod_into(tmp, uj[k][(t - i) % L], ft);
                }
                for (std::size_t h = 0; h < Mb.xpar.size(); ++h) {
                    if (t - Mb.xlag[h] < 0) continue;
                    value_jet(Mb.xZ[h], Mb.xK[h], t, tmp);
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
                for (int c = 0; c < mk; ++c) vrow[c] = seed[a](row, act[c] - 1);
                linear_jet(0.0, vrow.data(), dirs, Ej[a]);
            }
            for (int k = 0; k < K; ++k) {
                int a = M[k].eq;
                J5& Ea = Ej[a];
                const J5& F = fj[k][ts];
                for (int c = 0; c < mk; ++c) Ea.d1[c] += F.d1[c];
                Ea.M = F.M;
                Ea.linear = false;
                Ea.finish(dirs);
                for (int c = 0; c < mk; ++c) J_[k](row, act[c] - 1) = Ea.d1[c];
                for (int d = 0; d < 3; ++d) {
                    const double* r1 = Ea.row(1 << d);
                    const double* r2 = Ea.row(pairs[d]);
                    for (int c = 0; c < mk; ++c) {
                        Ph[k][d](row, act[c] - 1) = r1[c];
                        Ps[k][d](row, act[c] - 1) = r2[c];
                    }
                }
                const double* r3 = Ea.row(7);
                for (int c = 0; c < mk; ++c) Dm[k](row, act[c] - 1) = r3[c];
                const double* top = F.mat(7);
                double gk = g(row, k);
                for (std::size_t e = 0; e < mk2; ++e) Wl[e] += gk * top[e];
            }
            for (int a = 0; a < P; ++a) E[a] = &Ej[a];
            for (int l = 0; l < K; ++l) {
                Contract C;
                C.A = &Ar;
                C.row = row;
                C.lofs = (std::size_t) l;
                C.kmult = (std::size_t) K;
                score_jet(C, E, u_at(row, l), uj[l][ts], dirs);
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
