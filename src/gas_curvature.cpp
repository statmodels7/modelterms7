#include <Rcpp.h>
#include <fenv.h>
#include <RcppParallel.h>
#include <vector>
using namespace Rcpp;

// The second-order recursion of the submodel route, compiled: the exact
// curvature of a score-driven predictor in the term's own parameters,
// carried per group on that group's active coordinates. This is the stage
// piano_parallel.txt's voce 5 left blocked on the R recursion; what makes
// it portable is that the callbacks here are LOOKUPS (the derivatives are
// read at a predictor the caller has already fitted, the regime asymmetry)
// and the layer's blocks callback is arithmetic over vectors the caller
// already holds, so the loop touches no R API and the GROUPS can run over
// threads -- each group's W_l is accumulated locally and the merge into W
// runs on the main thread in group order, so no reduction is split and the
// result does not depend on the count, bit for bit.
//
// COVERAGE, stated rather than blurred: second order only (a direction
// routes to the R recursion, which stays the reference and the third-order
// implementation), and constant autoregressive charts per group
// (a developed partial autocorrelation keeps the R route). Every
// expression MIRRORS the R route's, matrix addition by matrix addition
// and product by product, so on a compiler that does not contract the two
// routes agree to the bit; the twin test compares at a tolerance for the
// platforms whose compiler does (the seg_block rule).
//
// Each element of `groups` carries, all built by .gas_curv_prep():
//   rows (1-based, time order), act (1-based columns among the caller's m),
//   wom/za (m_g x mk scattered chained rows and raw chart rows per
//   parameter), k2 vectors, db rows, hb (q constant mk x mk blocks),
//   f0, f0_u, f0_uu.

namespace {

// The same shape as gas_filter.cpp's, and for the reasons stated there:
// one compiled copy for both branches, and the calling thread's
// floating-point environment installed before the chunk.
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

struct CurvGrp {
    const int* rows;
    const int* act;
    int m, mk;
    const double* wom;                  // m x mk
    const double* zom;                  // m x mk (raw chart rows)
    const double* k2om;                 // m
    std::vector<const double*> wa, za, k2a;   // p of each
    std::vector<const double*> db;      // q, m x mk
    std::vector<const double*> hb;      // q, mk x mk (constant per group)
    double f0;
    const double* f0u;                  // mk
    const double* f0uu;                 // mk x mk
    double* Wl;                         // mk x mk output, zeroed
};

constexpr int kMinGroupsPar = 8;

} // namespace

// [[Rcpp::export]]
List gas_curvature_sub_cpp(NumericVector eta, List groups, int p, int q,
                           NumericVector om, NumericMatrix A,
                           NumericMatrix B, int ap,
                           NumericVector s_at, NumericVector c_at,
                           NumericVector g,
                           NumericMatrix Hc, NumericMatrix D3m,
                           List Vs, NumericMatrix seed,
                           int threads = 1) {
    int n = eta.size();
    int m_full = seed.ncol();
    int npar = Vs.size();               // the model's distribution parameters
    NumericMatrix D(n, m_full);

    int ng = groups.size();
    std::vector<CurvGrp> grps(ng);
    std::vector<List> keep(ng);
    std::vector<std::vector<double>> Wls(ng);
    for (int l = 0; l < ng; ++l) {
        keep[l] = as<List>(groups[l]);
        List grp = keep[l];
        CurvGrp& G = grps[l];
        IntegerVector rows = grp["rows"];
        IntegerVector act = grp["act"];
        G.rows = rows.begin();
        G.act = act.begin();
        G.m = rows.size();
        G.mk = act.size();
        NumericMatrix wom = grp["wom"]; G.wom = wom.begin();
        NumericMatrix zom = grp["zom"]; G.zom = zom.begin();
        NumericVector k2om = grp["k2om"]; G.k2om = k2om.begin();
        List wa = grp["wa"], za = grp["za"], k2a = grp["k2a"];
        for (int i = 0; i < p; ++i) {
            NumericMatrix w = wa[i]; G.wa.push_back(w.begin());
            NumericMatrix z = za[i]; G.za.push_back(z.begin());
            NumericVector k2 = k2a[i]; G.k2a.push_back(k2.begin());
        }
        List db = grp["db"], hb = grp["hb"];
        for (int j = 0; j < q; ++j) {
            NumericMatrix d = db[j]; G.db.push_back(d.begin());
            NumericMatrix h = hb[j]; G.hb.push_back(h.begin());
        }
        G.f0 = as<double>(grp["f0"]);
        NumericVector f0u = grp["f0u"]; G.f0u = f0u.begin();
        NumericMatrix f0uu = grp["f0uu"]; G.f0uu = f0uu.begin();
        Wls[l].assign((std::size_t) G.mk * G.mk, 0.0);
        G.Wl = Wls[l].data();
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
    const double* seed_p = seed.begin();
    double* D_p = D.begin();
    int An = A.nrow(), Bn = B.nrow();
    std::vector<const double*> Vs_p(npar);
    for (int a = 0; a < npar; ++a) {
        NumericMatrix V = Vs[a];
        Vs_p[a] = V.nrow() ? V.begin() : nullptr;
    }

    auto run_group = [&](std::size_t gi) {
        const CurvGrp& G = grps[gi];
        int m = G.m, mk = G.mk;
        std::size_t mk2 = (std::size_t) mk * mk;
        std::vector<double> f(m, 0.0), s(m, 0.0);
        std::vector<double> F_(m * (std::size_t) mk, 0.0);
        std::vector<double> Sd(m * (std::size_t) mk, 0.0);
        std::vector<double> Phi(m * mk2, 0.0), Sdd(m * mk2, 0.0);
        std::vector<double> Ft(mk), Dt(mk), cross(mk), vr((std::size_t) npar * mk);
        std::vector<double> Pt(mk2), M(mk2);

        for (int t = 0; t < m; ++t) {
            int r = G.rows[t];               // 1-based
            double ft = om_p[r - 1];
            for (int k = 0; k < mk; ++k) Ft[k] = G.wom[k * m + t];
            // Pt starts as the chart curvature of the level:
            // k2[t] * (z_j * z_k), the product first as in the R h_of
            for (int j = 0; j < mk; ++j) {
                double zj = G.zom[j * m + t];
                for (int k = 0; k < mk; ++k) {
                    Pt[(std::size_t) k * mk + j] =
                        G.k2om[t] * (zj * G.zom[k * m + t]);
                }
            }

            for (int i = 1; i <= p; ++i) {
                double a_r = A_p[(i - 1) * An + (r - 1)];
                const double* wa = G.wa[i - 1];
                const double* za = G.za[i - 1];
                double k2 = G.k2a[i - 1][t];
                double s_lag = (t - i >= 0) ? s[t - i] : 0.0;
                const double* Sd_l = (t - i >= 0) ? &Sd[(std::size_t)(t - i) * mk]
                                                  : nullptr;
                const double* Sdd_l = (t - i >= 0) ? &Sdd[(std::size_t)(t - i) * mk2]
                                                   : nullptr;
                ft += a_r * s_lag;
                for (int k = 0; k < mk; ++k) {
                    double sd = Sd_l ? Sd_l[k] : 0.0;
                    Ft[k] += a_r * sd;
                    Ft[k] += s_lag * wa[k * m + t];
                }
                for (int j = 0; j < mk; ++j) {
                    double au_j = wa[j * m + t];
                    double sd_j = Sd_l ? Sd_l[j] : 0.0;
                    double za_j = za[j * m + t];
                    for (int k = 0; k < mk; ++k) {
                        std::size_t jk = (std::size_t) k * mk + j;
                        double sd_k = Sd_l ? Sd_l[k] : 0.0;
                        double acc = Pt[jk];
                        acc += a_r * (Sdd_l ? Sdd_l[jk] : 0.0);
                        acc += au_j * sd_k;
                        acc += sd_j * wa[k * m + t];
                        acc += s_lag * (k2 * (za_j * za[k * m + t]));
                        Pt[jk] = acc;
                    }
                }
            }

            for (int j2 = 1; j2 <= q; ++j2) {
                double b_r = B_p[(j2 - 1) * Bn + (r - 1)];
                const double* dbj = G.db[j2 - 1];
                const double* hbj = G.hb[j2 - 1];
                double f_lag = (t - j2 >= 0) ? f[t - j2] : G.f0;
                const double* F_l = (t - j2 >= 0) ? &F_[(std::size_t)(t - j2) * mk]
                                                  : G.f0u;
                const double* Phi_l = (t - j2 >= 0) ? &Phi[(std::size_t)(t - j2) * mk2]
                                                    : G.f0uu;
                ft += b_r * f_lag;
                for (int k = 0; k < mk; ++k) {
                    Ft[k] += b_r * F_l[k];
                    Ft[k] += f_lag * dbj[k * m + t];
                }
                for (int j = 0; j < mk; ++j) {
                    double bu_j = dbj[j * m + t];
                    double Fl_j = F_l[j];
                    for (int k = 0; k < mk; ++k) {
                        std::size_t jk = (std::size_t) k * mk + j;
                        double acc = Pt[jk];
                        acc += b_r * Phi_l[jk];
                        acc += bu_j * F_l[k];
                        acc += Fl_j * dbj[k * m + t];
                        acc += f_lag * hbj[jk];
                        Pt[jk] = acc;
                    }
                }
            }

            f[t] = ft;
            for (int k = 0; k < mk; ++k) F_[(std::size_t) t * mk + k] = Ft[k];
            for (std::size_t jk = 0; jk < mk2; ++jk) {
                Phi[(std::size_t) t * mk2 + jk] = Pt[jk];
            }

            double e_t = eta_p[r - 1] + ft;
            (void) e_t;                     // the lookups do not read it
            for (int k = 0; k < mk; ++k) {
                Dt[k] = seed_p[(std::size_t)(G.act[k] - 1) * n + (r - 1)] + Ft[k];
                D_p[(std::size_t)(G.act[k] - 1) * n + (r - 1)] = Dt[k];
            }
            for (std::size_t jk = 0; jk < mk2; ++jk) {
                G.Wl[jk] += g_p[r - 1] * Pt[jk];
            }

            s[t] = s_p[r - 1];              // drives the later lags
            double cv = c_p[r - 1];
            // the blocks of the caller: cross and M from H, D3 and Vs, the
            // same loops in the same order as .structural_blocks()
            for (int k = 0; k < mk; ++k) cross[k] = 0.0;
            for (int qq = 0; qq < npar; ++qq) {
                if (qq == ap - 1) continue;
                double h = Hc_p[(std::size_t) qq * n + (r - 1)];
                const double* Vq = Vs_p[qq];
                for (int k = 0; k < mk; ++k) {
                    cross[k] += h * Vq[(std::size_t)(G.act[k] - 1) * n + (r - 1)];
                }
            }
            for (int a = 0; a < npar; ++a) {
                double* v = &vr[(std::size_t) a * mk];
                if (a == ap - 1) {
                    for (int k = 0; k < mk; ++k) v[k] = Dt[k];
                } else {
                    const double* Va = Vs_p[a];
                    for (int k = 0; k < mk; ++k) {
                        v[k] = Va[(std::size_t)(G.act[k] - 1) * n + (r - 1)];
                    }
                }
            }
            for (std::size_t jk = 0; jk < mk2; ++jk) M[jk] = 0.0;
            for (int a = 0; a < npar; ++a) {
                const double* va = &vr[(std::size_t) a * mk];
                for (int b = 0; b < npar; ++b) {
                    double d3 = D3_p[(std::size_t)(a * npar + b) * n + (r - 1)];
                    const double* vb = &vr[(std::size_t) b * mk];
                    for (int j = 0; j < mk; ++j) {
                        double vaj = va[j];
                        for (int k = 0; k < mk; ++k) {
                            M[(std::size_t) k * mk + j] += d3 * (vaj * vb[k]);
                        }
                    }
                }
            }
            for (int k = 0; k < mk; ++k) {
                Sd[(std::size_t) t * mk + k] = cv * Dt[k] + cross[k];
            }
            for (std::size_t jk = 0; jk < mk2; ++jk) {
                Sdd[(std::size_t) t * mk2 + jk] = cv * Pt[jk] + M[jk];
            }
        }
    };

    // the count is passed on rather than left to the process-level setting
    auto body = [&](std::size_t gi) { run_group(gi); };
    GroupWorker<decltype(body)> w(body);
    if (threads > 1 && ng >= kMinGroupsPar) {
        RcppParallel::parallelFor(0, (std::size_t) ng, w, 1, threads);
    } else {
        w(0, (std::size_t) ng);
    }

    // the merge into W runs HERE, on the main thread, in group order
    NumericMatrix W(m_full, m_full);
    double* W_p = W.begin();
    for (int l = 0; l < ng; ++l) {
        const CurvGrp& G = grps[l];
        for (int j = 0; j < G.mk; ++j) {
            for (int k = 0; k < G.mk; ++k) {
                W_p[(std::size_t)(G.act[k] - 1) * m_full + (G.act[j] - 1)] +=
                    G.Wl[(std::size_t) k * G.mk + j];
            }
        }
    }
    return List::create(_["jacobian"] = D, _["curvature"] = W);
}
