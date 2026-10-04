#include <Rcpp.h>
#include <fenv.h>
#include <RcppParallel.h>
#include <vector>
using namespace Rcpp;

// The second, third and fourth orders of the score-driven recursion,
// compiled, for both routes of term_curvature(), term_third() and
// term_fourth(). It mirrors .gas_curvature_sub() expression by expression:
// the chart of every value of the recursion (the level, the loadings, the
// partial autocorrelations) is read per observation through its raw design
// row z and its link's derivatives k2, k3, k4 -- h = k2 z z', the third
// derivative along v is k3 (z.v) z z', the fourth along v and w is
// k4 (z.v)(z.w) z z' -- which is exact for a scalar parameter too, its z
// being the unit vector at its coordinate. The autoregressive coefficients
// go through the Levinson-Durbin map, whose derivatives to the fourth order
// are written here as gas_levinson2/3/4 write them in R. The model's own
// pieces (cross, M, dcurv, N, Q, P, cppp) are read from the arrays the
// caller supplies, as .structural_blocks() reads them.
//
// The lags are held in ring buffers of max(p, q) + 1 slots, so the memory
// is constant in the length of a group. Groups run over threads with each
// group's W accumulated locally and merged on the main thread in group
// order, so no reduction is split and the result does not depend on the
// thread count, bit for bit.

namespace {

#if defined(__GNUC__) || defined(__clang__)
#define MT7_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define MT7_NOINLINE __declspec(noinline)
#else
#define MT7_NOINLINE
#endif

template <typename Body>
struct GroupWorker : public RcppParallel::Worker {
    const Body& body;
    fenv_t env;
    explicit GroupWorker(const Body& b) : body(b) { fegetenv(&env); }
    MT7_NOINLINE void operator()(std::size_t begin, std::size_t end) {
        fesetenv(&env);
        for (std::size_t g = begin; g < end; ++g) body(g);
    }
};

constexpr int kMinGroupsPar = 8;

// one value of the recursion on a group: its per-row value, chained
// derivative row W, raw chart row Z (both m_g x mk, column-major) and the
// link's second to fourth derivatives per row
struct Val {
    const double* v = nullptr;
    const double* W = nullptr;
    const double* Z = nullptr;
    const double* k2 = nullptr;
    const double* k3 = nullptr;
    const double* k4 = nullptr;
};

struct Grp {
    const int* rows;
    const int* act;
    int m, mk;
    const double* vks;                 // mk x nd
    Val om;
    std::vector<Val> al;               // p
    std::vector<Val> pa;               // q, read only where b varies
    std::vector<const double*> db;     // q, m x mk chained rows of b
    bool varying_b;
    std::vector<const double*> hb;     // q, mk x mk, constant b only
    std::vector<const double*> tb;     // nd * q, mk x mk, constant b only
    std::vector<const double*> qb;     // q, mk x mk, constant b only
    double f0;
    const double* f0u;                 // mk
    const double* f0uu;                // mk x mk
    std::vector<const double*> f0_3;   // nd, mk x mk
    std::vector<double> df0;           // nd
    const double* dphi0;               // mk x nd
    const double* f0_4;                // mk x mk
};

// ---- small dense helpers on mk x mk column-major blocks ---------------------
inline void zero(std::vector<double>& x) { std::fill(x.begin(), x.end(), 0.0); }

// h += c * a b'
inline void add_outer(double* h, double c, const double* a, const double* b,
                      int mk) {
    for (int k = 0; k < mk; ++k) {
        double bk = b[k];
        for (int j = 0; j < mk; ++j) h[(std::size_t) k * mk + j] += c * (a[j] * bk);
    }
}
// h += a b' + b a'
inline void add_sym_outer(double* h, const double* a, const double* b, int mk) {
    for (int k = 0; k < mk; ++k) {
        for (int j = 0; j < mk; ++j) {
            h[(std::size_t) k * mk + j] += a[j] * b[k] + b[j] * a[k];
        }
    }
}
inline void matvec(const double* H, const double* v, double* out, int mk) {
    for (int j = 0; j < mk; ++j) out[j] = 0.0;
    for (int k = 0; k < mk; ++k) {
        double vk = v[k];
        for (int j = 0; j < mk; ++j) out[j] += H[(std::size_t) k * mk + j] * vk;
    }
}
inline double dot(const double* a, const double* b, int mk) {
    double s = 0.0;
    for (int k = 0; k < mk; ++k) s += a[k] * b[k];
    return s;
}

// the chart's per-row pieces of one value: h = k2 z z', the third along d,
// the fourth along both directions
inline void chart_h(const Val& V, int t, int m, int mk, double c, double* h) {
    double k2 = V.k2[t];
    if (k2 == 0.0 || c == 0.0) return;
    for (int k = 0; k < mk; ++k) {
        double zk = V.Z[(std::size_t) k * m + t];
        if (zk == 0.0) continue;
        for (int j = 0; j < mk; ++j) {
            double zj = V.Z[(std::size_t) j * m + t];
            h[(std::size_t) k * mk + j] += c * (k2 * (zj * zk));
        }
    }
}
inline double zdot(const Val& V, int t, int m, int mk, const double* v) {
    double s = 0.0;
    for (int k = 0; k < mk; ++k) s += V.Z[(std::size_t) k * m + t] * v[k];
    return s;
}
inline void chart_t(const Val& V, int t, int m, int mk, const double* vd,
                    double c, double* h) {
    double k3 = V.k3[t];
    if (k3 == 0.0 || c == 0.0) return;
    double s = k3 * zdot(V, t, m, mk, vd);
    for (int k = 0; k < mk; ++k) {
        double zk = V.Z[(std::size_t) k * m + t];
        if (zk == 0.0) continue;
        for (int j = 0; j < mk; ++j) {
            double zj = V.Z[(std::size_t) j * m + t];
            h[(std::size_t) k * mk + j] += c * (s * (zj * zk));
        }
    }
}
inline void chart_q(const Val& V, int t, int m, int mk, const double* v1,
                    const double* v2, double c, double* h) {
    double k4 = V.k4[t];
    if (k4 == 0.0 || c == 0.0) return;
    double s = k4 * zdot(V, t, m, mk, v1) * zdot(V, t, m, mk, v2);
    for (int k = 0; k < mk; ++k) {
        double zk = V.Z[(std::size_t) k * m + t];
        if (zk == 0.0) continue;
        for (int j = 0; j < mk; ++j) {
            double zj = V.Z[(std::size_t) j * m + t];
            h[(std::size_t) k * mk + j] += c * (s * (zj * zk));
        }
    }
}
inline void row_of(const Val& V, int t, int m, int mk, double* out) {
    for (int k = 0; k < mk; ++k) out[k] = V.W[(std::size_t) k * m + t];
}

// ---- the Levinson-Durbin map's derivatives, as gas_levinson2/3/4 ------------
// J is q x q (row j = coefficient, column k = pacf), H[j], T[j], Q[j] are
// q x q, all column-major
struct Lev {
    int q;
    std::vector<double> J, H, T1, T2, Qd;   // H, T1, T2, Qd: q blocks of q x q
};

// the second-order map with, optionally, the third along w1 and w2 and the
// fourth along (w1, w2)
inline void levinson(const double* rho, int q, const double* w1,
                     const double* w2, int order, Lev& L) {
    L.q = q;
    std::size_t qq = (std::size_t) q * q;
    std::vector<double> phi, jac, hes, t1, t2, qd;
    for (int k = 1; k <= q; ++k) {
        std::vector<double> nphi(k), njac((std::size_t) k * q, 0.0),
            nhes(k * qq, 0.0), nt1(k * qq, 0.0), nt2(k * qq, 0.0),
            nqd(k * qq, 0.0);
        double rk = rho[k - 1];
        nphi[k - 1] = rk;
        // njac is k x q column-major: (i, c) at c * k + i
        njac[(std::size_t)(k - 1) * k + (k - 1)] = 1.0;
        int km1 = k - 1;
        for (int i = 0; i < km1; ++i) {
            int r = km1 - 1 - i;               // rev index, 0-based
            nphi[i] = phi[i] - rk * phi[r];
            for (int c = 0; c < q; ++c) {
                njac[(std::size_t) c * k + i] =
                    jac[(std::size_t) c * km1 + i] - rk * jac[(std::size_t) c * km1 + r];
            }
            njac[(std::size_t)(k - 1) * k + i] -= phi[r];
            // the hessian of coefficient i
            double* h = &nhes[i * qq];
            const double* hi = &hes[i * qq];
            const double* hr = &hes[r * qq];
            for (std::size_t e = 0; e < qq; ++e) h[e] = hi[e] - rk * hr[e];
            for (int c = 0; c < q; ++c) {
                double jr = jac[(std::size_t) c * km1 + r];
                h[(std::size_t) c * q + (k - 1)] -= jr;     // row k
                h[(std::size_t)(k - 1) * q + c] -= jr;      // column k
            }
            if (order >= 3) {
                // third along w1 (and w2)
                for (int which = 0; which < (order >= 4 || w2 ? 2 : 1); ++which) {
                    const double* w = which == 0 ? w1 : w2;
                    if (!w) continue;
                    std::vector<double>& nt = which == 0 ? nt1 : nt2;
                    const std::vector<double>& tt = which == 0 ? t1 : t2;
                    double* T = &nt[i * qq];
                    std::vector<double> hw(q, 0.0);
                    for (int c = 0; c < q; ++c) {
                        for (int b = 0; b < q; ++b) hw[b] += hr[(std::size_t) c * q + b] * w[c];
                    }
                    for (std::size_t e = 0; e < qq; ++e) {
                        T[e] = tt[i * qq + e] - rk * tt[r * qq + e] - w[k - 1] * hr[e];
                    }
                    for (int c = 0; c < q; ++c) {
                        T[(std::size_t) c * q + (k - 1)] -= hw[c];
                        T[(std::size_t)(k - 1) * q + c] -= hw[c];
                    }
                }
            }
            if (order >= 4) {
                double* Qm = &nqd[i * qq];
                const double* t1r = &t1[r * qq];
                const double* t2r = &t2[r * qq];
                std::vector<double> tvw(q, 0.0);
                for (int c = 0; c < q; ++c) {
                    for (int b = 0; b < q; ++b) tvw[b] += t1r[(std::size_t) c * q + b] * w2[c];
                }
                for (std::size_t e = 0; e < qq; ++e) {
                    Qm[e] = qd[i * qq + e] - rk * qd[r * qq + e] -
                        w1[k - 1] * t2r[e] - w2[k - 1] * t1r[e];
                }
                for (int c = 0; c < q; ++c) {
                    Qm[(std::size_t) c * q + (k - 1)] -= tvw[c];
                    Qm[(std::size_t)(k - 1) * q + c] -= tvw[c];
                }
            }
        }
        phi.swap(nphi); jac.swap(njac); hes.swap(nhes);
        t1.swap(nt1); t2.swap(nt2); qd.swap(nqd);
    }
    L.J = jac;           // q x q, (j, c) at c * q + j
    L.H = hes;
    L.T1 = t1;
    L.T2 = t2;
    L.Qd = qd;
}

} // namespace

// [[Rcpp::export]]
List gas_curvature_gen_cpp(NumericVector eta, List groups, int p, int q,
                           int nd, NumericVector om, NumericMatrix A,
                           NumericMatrix B, int ap, NumericVector s_at,
                           NumericVector c_at, NumericVector g,
                           NumericMatrix Hc, NumericMatrix D3m,
                           NumericMatrix D4m, NumericMatrix D5m,
                           List Vs, NumericMatrix seed, int threads = 1) {
    int n = eta.size();
    int m_full = seed.ncol();
    int npar = Vs.size();
    bool third = nd >= 1, fourth = nd >= 2;
    NumericMatrix D(n, m_full);
    std::vector<NumericMatrix> dP;
    for (int d = 0; d < nd; ++d) dP.push_back(NumericMatrix(n, m_full));
    NumericMatrix dPS(fourth ? n : 0, fourth ? m_full : 0);

    int ng = groups.size();
    std::vector<Grp> grps(ng);
    std::vector<List> keep(ng);
    std::vector<std::vector<double>> Wls(ng);
    auto read_val = [&](List x) {
        Val V;
        NumericVector v = x["v"]; V.v = v.begin();
        NumericMatrix W = x["W"]; V.W = W.begin();
        NumericMatrix Z = x["Z"]; V.Z = Z.begin();
        NumericVector k2 = x["k2"]; V.k2 = k2.begin();
        if (third) { NumericVector k3 = x["k3"]; V.k3 = k3.begin(); }
        if (fourth) { NumericVector k4 = x["k4"]; V.k4 = k4.begin(); }
        return V;
    };
    for (int l = 0; l < ng; ++l) {
        keep[l] = as<List>(groups[l]);
        List grp = keep[l];
        Grp& G = grps[l];
        IntegerVector rows = grp["rows"]; G.rows = rows.begin(); G.m = rows.size();
        IntegerVector act = grp["act"]; G.act = act.begin(); G.mk = act.size();
        if (third) { NumericMatrix vks = grp["vks"]; G.vks = vks.begin(); }
        G.om = read_val(as<List>(grp["om"]));
        List al = grp["al"];
        for (int i = 0; i < p; ++i) G.al.push_back(read_val(as<List>(al[i])));
        G.varying_b = as<bool>(grp["varying_b"]);
        if (q > 0) {
            List db = grp["db"];
            for (int j = 0; j < q; ++j) { NumericMatrix x = db[j]; G.db.push_back(x.begin()); }
            if (G.varying_b) {
                List pa = grp["pa"];
                for (int j = 0; j < q; ++j) G.pa.push_back(read_val(as<List>(pa[j])));
            } else {
                List hb = grp["hb"];
                for (int j = 0; j < q; ++j) { NumericMatrix x = hb[j]; G.hb.push_back(x.begin()); }
                if (third) {
                    List tb = grp["tb"];
                    for (int d = 0; d < nd; ++d) {
                        List tbd = tb[d];
                        for (int j = 0; j < q; ++j) { NumericMatrix x = tbd[j]; G.tb.push_back(x.begin()); }
                    }
                }
                if (fourth) {
                    List qb = grp["qb"];
                    for (int j = 0; j < q; ++j) { NumericMatrix x = qb[j]; G.qb.push_back(x.begin()); }
                }
            }
        }
        G.f0 = as<double>(grp["f0"]);
        NumericVector f0u = grp["f0u"]; G.f0u = f0u.begin();
        NumericMatrix f0uu = grp["f0uu"]; G.f0uu = f0uu.begin();
        if (third) {
            List f3 = grp["f0_3"];
            for (int d = 0; d < nd; ++d) { NumericMatrix x = f3[d]; G.f0_3.push_back(x.begin()); }
            NumericVector df0 = grp["df0"];
            G.df0.assign(df0.begin(), df0.end());
            NumericMatrix dphi0 = grp["dphi0"]; G.dphi0 = dphi0.begin();
        }
        if (fourth) { NumericMatrix f04 = grp["f0_4"]; G.f0_4 = f04.begin(); }
        Wls[l].assign((std::size_t) G.mk * G.mk, 0.0);
    }

    const double* eta_p = eta.begin();
    const double* om_p = om.begin();
    const double* A_p = A.begin();
    const double* B_p = B.begin();
    const double* s_p = s_at.begin();
    const double* c_p = c_at.begin();
    const double* g_p = g.begin();
    const double* Hc_p = Hc.begin();
    const double* D3_p = D3m.begin();
    const double* D4_p = D4m.nrow() ? D4m.begin() : nullptr;
    const double* D5_p = D5m.nrow() ? D5m.begin() : nullptr;
    const double* seed_p = seed.begin();
    double* D_p = D.begin();
    std::vector<double*> dP_p(nd);
    for (int d = 0; d < nd; ++d) dP_p[d] = dP[d].begin();
    double* dPS_p = fourth ? dPS.begin() : nullptr;
    int An = A.nrow(), Bn = B.nrow();
    std::vector<const double*> Vs_p(npar);
    for (int a = 0; a < npar; ++a) {
        NumericMatrix V = Vs[a];
        Vs_p[a] = V.nrow() ? V.begin() : nullptr;
    }
    std::size_t np2 = (std::size_t) npar * npar, np3 = np2 * npar;

    auto run_group = [&](std::size_t gi) {
        const Grp& G = grps[gi];
        int m = G.m, mk = G.mk;
        std::size_t mk2 = (std::size_t) mk * mk;
        int L = std::max(p, q) + 1;          // ring buffer slots
        auto slot = [&](int t) { return (std::size_t)(t % L); };
        std::vector<double> f(L, 0.0), s(L, 0.0);
        std::vector<double> F_(L * (std::size_t) mk, 0.0), Sd(L * (std::size_t) mk, 0.0);
        std::vector<double> Phi(L * mk2, 0.0), Sdd(L * mk2, 0.0);
        std::vector<double> Psi(third ? L * mk2 * nd : 0, 0.0),
            Sddd(third ? L * mk2 * nd : 0, 0.0);
        std::vector<double> Om(fourth ? L * mk2 : 0, 0.0), S4(fourth ? L * mk2 : 0, 0.0);
        std::vector<double> Ft(mk), Dt(mk), cross(mk), dcurv(mk), vr((std::size_t) npar * mk);
        std::vector<double> Pt(mk2), M(mk2), Qb(mk2), Pb(mk2);
        std::vector<double> Tt(third ? mk2 * nd : 0), Qt(fourth ? mk2 : 0);
        std::vector<double> Nd(third ? mk2 * nd : 0);
        std::vector<double> a_u(mk), tmp(mk), tmp2(mk), dvs((std::size_t) std::max(nd, 1) * npar);
        std::vector<double> Wl(mk2, 0.0);
        std::vector<double> zvec(mk, 0.0), zmat(mk2, 0.0);
        // per-row b pieces where b varies
        std::vector<double> hbv(q * mk2), tbv((std::size_t) q * mk2 * std::max(nd, 1)), qbv(q * mk2);
        std::vector<double> dbsum(mk);
        const double* v1 = third ? G.vks : nullptr;
        const double* v2 = fourth ? G.vks + mk : nullptr;

        // a jet: value, d1, d2, d3v, d3w, d4 and the contractions of the
        // lower orders, as .gas_jet()
        struct Jet {
            double v, a, b, d2vw;
            const double *d1, *d2, *d3v, *d3w, *d4;
            std::vector<double> d2v, d2w, d3vw;
        };
        auto make_jet = [&](double v0, const double* d1, const double* d2,
                            const double* d3v, const double* d3w,
                            const double* d4, Jet& J) {
            J.v = v0; J.d1 = d1; J.d2 = d2; J.d3v = d3v; J.d3w = d3w; J.d4 = d4;
            J.d2v.assign(mk, 0.0); J.d2w.assign(mk, 0.0); J.d3vw.assign(mk, 0.0);
            matvec(d2, v2, J.d2w.data(), mk);
            J.a = dot(d1, v1, mk);
            J.b = dot(d1, v2, mk);
            matvec(d2, v1, J.d2v.data(), mk);
            J.d2vw = dot(v1, J.d2w.data(), mk);
            matvec(d3v, v2, J.d3vw.data(), mk);
        };
        // .gas_prod4(A, B) added into out
        auto prod4 = [&](const Jet& A_, const Jet& B_, double* out) {
            for (std::size_t e = 0; e < mk2; ++e) {
                out[e] += A_.v * B_.d4[e] + A_.a * B_.d3w[e] + A_.b * B_.d3v[e];
            }
            add_sym_outer(out, A_.d1, B_.d3vw.data(), mk);
            for (std::size_t e = 0; e < mk2; ++e) out[e] += A_.d2vw * B_.d2[e];
            add_sym_outer(out, A_.d2v.data(), B_.d2w.data(), mk);
            add_sym_outer(out, A_.d2w.data(), B_.d2v.data(), mk);
            for (std::size_t e = 0; e < mk2; ++e) out[e] += A_.d2[e] * B_.d2vw;
            add_sym_outer(out, A_.d3vw.data(), B_.d1, mk);
            for (std::size_t e = 0; e < mk2; ++e) {
                out[e] += A_.d3v[e] * B_.b + A_.d3w[e] * B_.a + A_.d4[e] * B_.v;
            }
        };

        std::vector<double> rho(q), wdot1(q), wdot2(q), hvw(q);
        std::vector<std::vector<double>> wrows(q, std::vector<double>(mk)),
            hrows(q, std::vector<double>(mk2)), hv1(q, std::vector<double>(mk)),
            hv2(q, std::vector<double>(mk)), trow1(q, std::vector<double>(mk2)),
            trow2(q, std::vector<double>(mk2)), qrow(q, std::vector<double>(mk2)),
            tvw(q, std::vector<double>(mk));
        Lev Lv, Lw, Lm, L4;

        // the b pieces at one row where the partial autocorrelations vary,
        // as hb_of(), tb_of() and qb_of()
        auto b_pieces = [&](int t) {
            for (int k = 0; k < q; ++k) {
                const Val& V = G.pa[k];
                rho[k] = V.v[t];
                row_of(V, t, m, mk, wrows[k].data());
                zero(hrows[k]);
                chart_h(V, t, m, mk, 1.0, hrows[k].data());
                if (third) {
                    zero(trow1[k]);
                    chart_t(V, t, m, mk, v1, 1.0, trow1[k].data());
                    wdot1[k] = dot(wrows[k].data(), v1, mk);
                    matvec(hrows[k].data(), v1, hv1[k].data(), mk);
                }
                if (fourth) {
                    zero(trow2[k]);
                    chart_t(V, t, m, mk, v2, 1.0, trow2[k].data());
                    wdot2[k] = dot(wrows[k].data(), v2, mk);
                    matvec(hrows[k].data(), v2, hv2[k].data(), mk);
                    zero(qrow[k]);
                    chart_q(V, t, m, mk, v1, v2, 1.0, qrow[k].data());
                    hvw[k] = dot(v1, hv2[k].data(), mk);
                    matvec(trow1[k].data(), v2, tvw[k].data(), mk);
                }
            }
            std::size_t qq = (std::size_t) q * q;
            levinson(rho.data(), q, third ? wdot1.data() : nullptr,
                     fourth ? wdot2.data() : nullptr, fourth ? 4 : (third ? 3 : 2), Lv);
            if (fourth) levinson(rho.data(), q, hvw.data(), nullptr, 3, Lm);
            for (int j = 0; j < q; ++j) {
                const double* Hj = &Lv.H[j * qq];
                double* h = &hbv[j * mk2];
                for (std::size_t e = 0; e < mk2; ++e) h[e] = 0.0;
                for (int k = 0; k < q; ++k) {
                    for (int ll = 0; ll < q; ++ll) {
                        double c = Hj[(std::size_t) ll * q + k];
                        if (c != 0.0) add_outer(h, c, wrows[k].data(), wrows[ll].data(), mk);
                    }
                    double jc = Lv.J[(std::size_t) k * q + j];
                    if (jc != 0.0) {
                        for (std::size_t e = 0; e < mk2; ++e) h[e] += jc * hrows[k][e];
                    }
                }
                for (int d = 0; d < nd; ++d) {
                    const double* vd = d == 0 ? v1 : v2;
                    const std::vector<double>& wd = d == 0 ? wdot1 : wdot2;
                    const std::vector<std::vector<double>>& hvd = d == 0 ? hv1 : hv2;
                    const std::vector<std::vector<double>>& trd = d == 0 ? trow1 : trow2;
                    const double* T3 = d == 0 ? &Lv.T1[j * qq] : &Lv.T2[j * qq];
                    double* tbj = &tbv[((std::size_t) d * q + j) * mk2];
                    for (std::size_t e = 0; e < mk2; ++e) tbj[e] = 0.0;
                    std::vector<double> hw(q, 0.0);
                    for (int c = 0; c < q; ++c) {
                        for (int b = 0; b < q; ++b) hw[b] += Hj[(std::size_t) c * q + b] * wd[c];
                    }
                    for (int k = 0; k < q; ++k) {
                        for (int ll = 0; ll < q; ++ll) {
                            double c3 = T3[(std::size_t) ll * q + k];
                            if (c3 != 0.0) add_outer(tbj, c3, wrows[k].data(), wrows[ll].data(), mk);
                            double c2 = Hj[(std::size_t) ll * q + k];
                            if (c2 != 0.0) {
                                add_outer(tbj, c2, hvd[k].data(), wrows[ll].data(), mk);
                                add_outer(tbj, c2, wrows[k].data(), hvd[ll].data(), mk);
                            }
                        }
                        if (hw[k] != 0.0) {
                            for (std::size_t e = 0; e < mk2; ++e) tbj[e] += hw[k] * hrows[k][e];
                        }
                        double jc = Lv.J[(std::size_t) k * q + j];
                        if (jc != 0.0) {
                            for (std::size_t e = 0; e < mk2; ++e) tbj[e] += jc * trd[k][e];
                        }
                    }
                    (void) vd;
                }
                if (fourth) {
                    // the fourth order of b_j, as qb_of(): ld4 and ld3 along hvw
                    const double* Q4 = &Lv.Qd[j * qq];
                    const double* T3m = &Lm.T1[j * qq];
                    const double* T3v = &Lv.T1[j * qq];
                    const double* T3w = &Lv.T2[j * qq];
                    double* qbj = &qbv[j * mk2];
                    for (std::size_t e = 0; e < mk2; ++e) qbj[e] = 0.0;
                    std::vector<double> hv4(q, 0.0), hw4(q, 0.0), hm(q, 0.0), t3vw(q, 0.0);
                    for (int c = 0; c < q; ++c) {
                        for (int b = 0; b < q; ++b) {
                            double hjbc = Hj[(std::size_t) c * q + b];
                            hv4[b] += hjbc * wdot1[c];
                            hw4[b] += hjbc * wdot2[c];
                            hm[b] += hjbc * hvw[c];
                            t3vw[b] += T3v[(std::size_t) c * q + b] * wdot2[c];
                        }
                    }
                    for (int k = 0; k < q; ++k) {
                        for (int ll = 0; ll < q; ++ll) {
                            std::size_t kl = (std::size_t) ll * q + k;
                            double co = Q4[kl] + T3m[kl];
                            if (co != 0.0) add_outer(qbj, co, wrows[k].data(), wrows[ll].data(), mk);
                            if (T3v[kl] != 0.0) {
                                add_outer(qbj, T3v[kl], hv2[k].data(), wrows[ll].data(), mk);
                                add_outer(qbj, T3v[kl], wrows[k].data(), hv2[ll].data(), mk);
                            }
                            if (T3w[kl] != 0.0) {
                                add_outer(qbj, T3w[kl], hv1[k].data(), wrows[ll].data(), mk);
                                add_outer(qbj, T3w[kl], wrows[k].data(), hv1[ll].data(), mk);
                            }
                            double c2 = Hj[kl];
                            if (c2 != 0.0) {
                                add_outer(qbj, c2, tvw[k].data(), wrows[ll].data(), mk);
                                add_outer(qbj, c2, wrows[k].data(), tvw[ll].data(), mk);
                                add_outer(qbj, c2, hv1[k].data(), hv2[ll].data(), mk);
                                add_outer(qbj, c2, hv2[k].data(), hv1[ll].data(), mk);
                            }
                        }
                        if (t3vw[k] != 0.0) for (std::size_t e = 0; e < mk2; ++e) qbj[e] += t3vw[k] * hrows[k][e];
                        if (hm[k] != 0.0) for (std::size_t e = 0; e < mk2; ++e) qbj[e] += hm[k] * hrows[k][e];
                        if (hw4[k] != 0.0) for (std::size_t e = 0; e < mk2; ++e) qbj[e] += hw4[k] * trow1[k][e];
                        if (hv4[k] != 0.0) for (std::size_t e = 0; e < mk2; ++e) qbj[e] += hv4[k] * trow2[k][e];
                        double jc = Lv.J[(std::size_t) k * q + j];
                        if (jc != 0.0) for (std::size_t e = 0; e < mk2; ++e) qbj[e] += jc * qrow[k][e];
                    }
                }
            }
        };

        for (int t = 0; t < m; ++t) {
            int r = G.rows[t];
            std::size_t ts = slot(t);
            double ft = om_p[r - 1];
            row_of(G.om, t, m, mk, Ft.data());
            zero(Pt);
            chart_h(G.om, t, m, mk, 1.0, Pt.data());
            if (third) {
                for (int d = 0; d < nd; ++d) {
                    double* T = &Tt[(std::size_t) d * mk2];
                    for (std::size_t e = 0; e < mk2; ++e) T[e] = 0.0;
                    chart_t(G.om, t, m, mk, d == 0 ? v1 : v2, 1.0, T);
                }
            }
            if (fourth) {
                zero(Qt);
                chart_q(G.om, t, m, mk, v1, v2, 1.0, Qt.data());
            }

            for (int i = 1; i <= p; ++i) {
                int lag = t - i;
                bool has = lag >= 0;
                std::size_t ls = has ? slot(lag) : 0;
                double s_l = has ? s[ls] : 0.0;
                const double* Sd_l = has ? &Sd[ls * mk] : zvec.data();
                const double* Sdd_l = has ? &Sdd[ls * mk2] : zmat.data();
                const Val& V = G.al[i - 1];
                double av = V.v[t];
                row_of(V, t, m, mk, a_u.data());
                ft += av * s_l;
                for (int k = 0; k < mk; ++k) Ft[k] += av * Sd_l[k] + s_l * a_u[k];
                for (std::size_t e = 0; e < mk2; ++e) Pt[e] += av * Sdd_l[e];
                add_sym_outer(Pt.data(), a_u.data(), Sd_l, mk);
                chart_h(V, t, m, mk, s_l, Pt.data());
                if (third) {
                    for (int d = 0; d < nd; ++d) {
                        const double* vd = d == 0 ? v1 : v2;
                        double* T = &Tt[(std::size_t) d * mk2];
                        const double* Sddd_l = has ? &Sddd[(ls * nd + d) * mk2] : zmat.data();
                        double dSd_l = has ? dot(Sd_l, vd, mk) : 0.0;
                        if (has) matvec(Sdd_l, vd, tmp.data(), mk); else std::fill(tmp.begin(), tmp.end(), 0.0);
                        double da = dot(a_u.data(), vd, mk);
                        // da2 = a_uu vd, a_uu = k2 z z'
                        std::fill(tmp2.begin(), tmp2.end(), 0.0);
                        {
                            double k2 = V.k2[t];
                            if (k2 != 0.0) {
                                double zv = zdot(V, t, m, mk, vd);
                                for (int j = 0; j < mk; ++j) {
                                    tmp2[j] = k2 * V.Z[(std::size_t) j * m + t] * zv;
                                }
                            }
                        }
                        for (std::size_t e = 0; e < mk2; ++e) T[e] += da * Sdd_l[e] + av * Sddd_l[e];
                        add_sym_outer(T, tmp2.data(), Sd_l, mk);
                        add_sym_outer(T, a_u.data(), tmp.data(), mk);
                        chart_h(V, t, m, mk, dSd_l, T);
                        chart_t(V, t, m, mk, vd, s_l, T);
                    }
                }
                if (fourth) {
                    std::vector<double> au2(mk2, 0.0), at1(mk2, 0.0), at2(mk2, 0.0), aq(mk2, 0.0);
                    chart_h(V, t, m, mk, 1.0, au2.data());
                    chart_t(V, t, m, mk, v1, 1.0, at1.data());
                    chart_t(V, t, m, mk, v2, 1.0, at2.data());
                    chart_q(V, t, m, mk, v1, v2, 1.0, aq.data());
                    Jet Aj, Bj;
                    make_jet(av, a_u.data(), au2.data(), at1.data(), at2.data(), aq.data(), Aj);
                    if (has) {
                        make_jet(s_l, Sd_l, Sdd_l, &Sddd[(ls * nd + 0) * mk2],
                                 &Sddd[(ls * nd + 1) * mk2], &S4[ls * mk2], Bj);
                    } else {
                        make_jet(0.0, zvec.data(), zmat.data(), zmat.data(), zmat.data(), zmat.data(), Bj);
                    }
                    prod4(Aj, Bj, Qt.data());
                }
            }

            if (q > 0) {
                if (G.varying_b) b_pieces(t);
                for (int j = 1; j <= q; ++j) {
                    int lag = t - j;
                    bool has = lag >= 0;
                    std::size_t ls = has ? slot(lag) : 0;
                    double f_l = has ? f[ls] : G.f0;
                    const double* F_l = has ? &F_[ls * mk] : G.f0u;
                    const double* Phi_l = has ? &Phi[ls * mk2] : G.f0uu;
                    const double* bu = G.db[j - 1];
                    for (int k = 0; k < mk; ++k) tmp[k] = bu[(std::size_t) k * m + t];
                    double bv = B_p[(std::size_t)(j - 1) * Bn + (r - 1)];
                    const double* hb = G.varying_b ? &hbv[(j - 1) * mk2] : G.hb[j - 1];
                    ft += bv * f_l;
                    for (int k = 0; k < mk; ++k) Ft[k] += bv * F_l[k] + f_l * tmp[k];
                    for (std::size_t e = 0; e < mk2; ++e) Pt[e] += bv * Phi_l[e];
                    add_sym_outer(Pt.data(), tmp.data(), F_l, mk);
                    for (std::size_t e = 0; e < mk2; ++e) Pt[e] += f_l * hb[e];
                    if (third) {
                        for (int d = 0; d < nd; ++d) {
                            const double* vd = d == 0 ? v1 : v2;
                            double* T = &Tt[(std::size_t) d * mk2];
                            const double* Psi_l = has ? &Psi[(ls * nd + d) * mk2] : G.f0_3[d];
                            double dF_l = has ? dot(F_l, vd, mk) : G.df0[d];
                            std::vector<double> dPhi_l(mk);
                            if (has) matvec(Phi_l, vd, dPhi_l.data(), mk);
                            else for (int k = 0; k < mk; ++k) dPhi_l[k] = G.dphi0[(std::size_t) d * mk + k];
                            double db = dot(tmp.data(), vd, mk);
                            matvec(hb, vd, tmp2.data(), mk);
                            const double* tb = G.varying_b ? &tbv[((std::size_t) d * q + (j - 1)) * mk2]
                                                           : G.tb[(std::size_t) d * q + (j - 1)];
                            for (std::size_t e = 0; e < mk2; ++e) T[e] += db * Phi_l[e] + bv * Psi_l[e];
                            add_sym_outer(T, tmp2.data(), F_l, mk);
                            add_sym_outer(T, tmp.data(), dPhi_l.data(), mk);
                            for (std::size_t e = 0; e < mk2; ++e) T[e] += dF_l * hb[e] + f_l * tb[e];
                        }
                    }
                    if (fourth) {
                        const double* tb1 = G.varying_b ? &tbv[((std::size_t) 0 * q + (j - 1)) * mk2] : G.tb[(j - 1)];
                        const double* tb2 = G.varying_b ? &tbv[((std::size_t) 1 * q + (j - 1)) * mk2] : G.tb[(std::size_t) q + (j - 1)];
                        const double* qb = G.varying_b ? &qbv[(j - 1) * mk2] : G.qb[j - 1];
                        Jet Bj, Fj;
                        make_jet(bv, tmp.data(), hb, tb1, tb2, qb, Bj);
                        if (has) {
                            make_jet(f_l, F_l, Phi_l, &Psi[(ls * nd + 0) * mk2], &Psi[(ls * nd + 1) * mk2],
                                     &Om[ls * mk2], Fj);
                        } else {
                            make_jet(G.f0, G.f0u, G.f0uu, G.f0_3[0], G.f0_3[1], G.f0_4, Fj);
                        }
                        prod4(Bj, Fj, Qt.data());
                    }
                }
            }

            f[ts] = ft;
            for (int k = 0; k < mk; ++k) F_[ts * mk + k] = Ft[k];
            std::copy(Pt.begin(), Pt.end(), Phi.begin() + ts * mk2);
            if (third) {
                for (int d = 0; d < nd; ++d) {
                    std::copy(Tt.begin() + (std::size_t) d * mk2, Tt.begin() + (std::size_t)(d + 1) * mk2,
                              Psi.begin() + (ts * nd + d) * mk2);
                }
            }
            if (fourth) std::copy(Qt.begin(), Qt.end(), Om.begin() + ts * mk2);

            for (int k = 0; k < mk; ++k) {
                Dt[k] = seed_p[(std::size_t)(G.act[k] - 1) * n + (r - 1)] + Ft[k];
                D_p[(std::size_t)(G.act[k] - 1) * n + (r - 1)] = Dt[k];
            }
            const double* contrib = fourth ? Qt.data() : (third ? Tt.data() : Pt.data());
            for (std::size_t e = 0; e < mk2; ++e) Wl[e] += g_p[r - 1] * contrib[e];

            s[ts] = s_p[r - 1];
            double cv = c_p[r - 1];
            // the model's pieces, as .structural_blocks()
            for (int k = 0; k < mk; ++k) cross[k] = 0.0;
            for (int qq = 0; qq < npar; ++qq) {
                if (qq == ap - 1) continue;
                double h = Hc_p[(std::size_t) qq * n + (r - 1)];
                const double* Vq = Vs_p[qq];
                for (int k = 0; k < mk; ++k) cross[k] += h * Vq[(std::size_t)(G.act[k] - 1) * n + (r - 1)];
            }
            for (int a = 0; a < npar; ++a) {
                double* v = &vr[(std::size_t) a * mk];
                if (a == ap - 1) {
                    for (int k = 0; k < mk; ++k) v[k] = Dt[k];
                } else {
                    const double* Va = Vs_p[a];
                    for (int k = 0; k < mk; ++k) v[k] = Va[(std::size_t)(G.act[k] - 1) * n + (r - 1)];
                }
            }
            zero(M);
            for (int a = 0; a < npar; ++a) {
                for (int b = 0; b < npar; ++b) {
                    double d3 = D3_p[(std::size_t)(a * npar + b) * n + (r - 1)];
                    add_outer(M.data(), d3, &vr[(std::size_t) a * mk], &vr[(std::size_t) b * mk], mk);
                }
            }
            double* Sd_t = &Sd[ts * mk];
            double* Sdd_t = &Sdd[ts * mk2];
            for (int k = 0; k < mk; ++k) Sd_t[k] = cv * Dt[k] + cross[k];
            for (std::size_t e = 0; e < mk2; ++e) Sdd_t[e] = cv * Pt[e] + M[e];

            if (third) {
                for (int d = 0; d < nd; ++d) {
                    const double* vd = d == 0 ? v1 : v2;
                    for (int a = 0; a < npar; ++a) {
                        dvs[(std::size_t) d * npar + a] = dot(&vr[(std::size_t) a * mk], vd, mk);
                    }
                }
                for (int k = 0; k < mk; ++k) dcurv[k] = 0.0;
                for (int a = 0; a < npar; ++a) {
                    double c3 = D3_p[(std::size_t)((ap - 1) * npar + a) * n + (r - 1)];
                    for (int k = 0; k < mk; ++k) dcurv[k] += c3 * vr[(std::size_t) a * mk + k];
                }
                for (int d = 0; d < nd; ++d) {
                    double* N = &Nd[(std::size_t) d * mk2];
                    for (std::size_t e = 0; e < mk2; ++e) N[e] = 0.0;
                    for (int a = 0; a < npar; ++a) {
                        for (int b = 0; b < npar; ++b) {
                            double co = 0.0;
                            for (int c = 0; c < npar; ++c) {
                                co += D4_p[((std::size_t) a * np2 + (std::size_t) b * npar + c) * n + (r - 1)] *
                                    dvs[(std::size_t) d * npar + c];
                            }
                            if (co != 0.0) add_outer(N, co, &vr[(std::size_t) a * mk], &vr[(std::size_t) b * mk], mk);
                        }
                    }
                }
                for (int d = 0; d < nd; ++d) {
                    const double* vd = d == 0 ? v1 : v2;
                    matvec(Pt.data(), vd, tmp.data(), mk);
                    for (int k = 0; k < mk; ++k) dP_p[d][(std::size_t)(G.act[k] - 1) * n + (r - 1)] = tmp[k];
                    double dcv = dot(dcurv.data(), vd, mk);
                    double* S3 = &Sddd[(ts * nd + d) * mk2];
                    const double* T = &Tt[(std::size_t) d * mk2];
                    const double* N = &Nd[(std::size_t) d * mk2];
                    for (std::size_t e = 0; e < mk2; ++e) S3[e] = dcv * Pt[e] + cv * T[e] + N[e];
                    add_sym_outer(S3, tmp.data(), dcurv.data(), mk);
                }
            }
            if (fourth) {
                zero(Qb); zero(Pb);
                for (int a = 0; a < npar; ++a) {
                    for (int b = 0; b < npar; ++b) {
                        double qc = D4_p[((std::size_t)(ap - 1) * np2 + (std::size_t) a * npar + b) * n + (r - 1)];
                        if (qc != 0.0) add_outer(Qb.data(), qc, &vr[(std::size_t) a * mk], &vr[(std::size_t) b * mk], mk);
                        double pc = 0.0;
                        for (int c = 0; c < npar; ++c) {
                            for (int e2 = 0; e2 < npar; ++e2) {
                                pc += D5_p[((std::size_t) a * np3 + (std::size_t) b * np2 + (std::size_t) c * npar + e2) * n + (r - 1)] *
                                    dvs[c] * dvs[(std::size_t) npar + e2];
                            }
                        }
                        if (pc != 0.0) add_outer(Pb.data(), pc, &vr[(std::size_t) a * mk], &vr[(std::size_t) b * mk], mk);
                    }
                }
                double cppp = D3_p[(std::size_t)((ap - 1) * npar + (ap - 1)) * n + (r - 1)];
                std::vector<double> dPhi_v(mk), dPhi_w(mk), psi_vw(mk), dcv(mk), dcw(mk), rw(mk);
                matvec(Pt.data(), v1, dPhi_v.data(), mk);
                matvec(Pt.data(), v2, dPhi_w.data(), mk);
                double phi_vw = dot(v1, dPhi_w.data(), mk);
                matvec(&Tt[0], v2, psi_vw.data(), mk);
                for (int k = 0; k < mk; ++k) dPS_p[(std::size_t)(G.act[k] - 1) * n + (r - 1)] = psi_vw[k];
                matvec(Qb.data(), v1, dcv.data(), mk);
                matvec(Qb.data(), v2, dcw.data(), mk);
                double dc_v = dot(dcurv.data(), v1, mk);
                double dc_w = dot(dcurv.data(), v2, mk);
                double d2c = dot(v1, dcw.data(), mk) + cppp * phi_vw;
                for (int k = 0; k < mk; ++k) rw[k] = dcw[k] + cppp * dPhi_w[k];
                double* S4t = &S4[ts * mk2];
                const double* T1 = &Tt[0];
                const double* T2 = &Tt[mk2];
                for (std::size_t e = 0; e < mk2; ++e) {
                    S4t[e] = d2c * Pt[e] + dc_v * T2[e] + dc_w * T1[e] + cv * Qt[e] +
                        Pb[e] + phi_vw * Qb[e];
                }
                add_sym_outer(S4t, dPhi_w.data(), dcv.data(), mk);
                add_sym_outer(S4t, psi_vw.data(), dcurv.data(), mk);
                add_sym_outer(S4t, dPhi_v.data(), rw.data(), mk);
            }
        }
        std::copy(Wl.begin(), Wl.end(), Wls[gi].begin());
    };

    auto body = [&](std::size_t gi) { run_group(gi); };
    GroupWorker<decltype(body)> w(body);
    if (threads > 1 && ng >= kMinGroupsPar) {
        RcppParallel::parallelFor(0, (std::size_t) ng, w, 1, threads);
    } else {
        w(0, (std::size_t) ng);
    }

    NumericMatrix W(m_full, m_full);
    double* W_p = W.begin();
    for (int l = 0; l < ng; ++l) {
        const Grp& G = grps[l];
        for (int j = 0; j < G.mk; ++j) {
            for (int k = 0; k < G.mk; ++k) {
                W_p[(std::size_t)(G.act[k] - 1) * m_full + (G.act[j] - 1)] +=
                    Wls[l][(std::size_t) k * G.mk + j];
            }
        }
    }
    List dPl(nd);
    for (int d = 0; d < nd; ++d) dPl[d] = dP[d];
    return List::create(_["jacobian"] = D, _["curvature"] = W, _["dphi"] = dPl,
                        _["dpsi"] = dPS);
}
