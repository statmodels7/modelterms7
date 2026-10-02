#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

// The quadrature of a marginal seg() or jseg() term: the node set of one
// group, the shift each node adds to the group's observations, and the
// forward accumulation of the log-likelihood and its Jacobian. These are
// the bodies of .marg_nodes_r(), .marg_seg_shift_r() and the loop of
// .marg_seg_loglik_r() in R/marginal.R, written out; the family's density
// and score stay in R, as two vectorized calls per group made before the
// loop starts. The derivation of every formula is with those R twins, which
// the tests hold this file to.

// seq(from, to, length.out = n) as R's seq.default returns it: the two
// ends exact and the inside at from + i * by
static std::vector<double> seq_lo(double from, double to, int n) {
  std::vector<double> out(n);
  if (n == 1) {
    out[0] = from;
    return out;
  }
  const double by = (to - from) / (n - 1);
  out[0] = from;
  for (int i = 1; i < n - 1; ++i) out[i] = from + i * by;
  out[n - 1] = to;
  return out;
}

// [[Rcpp::export]]
List marg_seg_nodes_cpp(const NumericVector& xs, double m, double tau,
                        const NumericVector& gk_nodes,
                        const NumericVector& gk_wk, double mr) {
  const int nx = xs.size(), nn = gk_nodes.size();
  const double x1 = xs[0], xn = xs[nx - 1];
  const double eps = std::numeric_limits<double>::epsilon();
  std::vector<double> lo, hi;
  std::vector<int> edge;
  // the lower region
  const double a = std::min(m, x1) - 8.5 * tau;
  const double wid = x1 - a;
  // the panel counts are clamped in double before the cast, so a tiny tau
  // cannot overflow an int
  const int ns = (int) std::min(8.0, std::max(4.0, std::ceil(wid / (2.5 * tau))));
  std::vector<double> ed = seq_lo(a, x1, ns + 1);
  for (int s = 0; s < ns; ++s) {
    lo.push_back(ed[s]);
    hi.push_back(ed[s + 1]);
    edge.push_back(1);
  }
  // the interior intervals
  for (int j = 0; j < nx - 1; ++j) {
    const double wj = xs[j + 1] - xs[j];
    if (wj <= 0) continue;
    const int nsj = (int) std::min(6.0, std::max(1.0, std::ceil(wj / (2.5 * tau))));
    std::vector<double> edj = seq_lo(xs[j], xs[j + 1], nsj + 1);
    for (int s = 0; s < nsj; ++s) {
      lo.push_back(edj[s]);
      hi.push_back(edj[s + 1]);
      edge.push_back(0);
    }
  }
  const int np = lo.size();
  const int C = np * nn + 1;
  NumericVector p(C), lw(C), z(C), glw_m(C), glw_t(C), dpsi_m(C), dpsi_t(C),
    alw_mm(C), alw_mt(C), alw_tt(C);
  const double denom = std::max(x1 - a, eps);
  const double da_m = m < x1 ? 1.0 : 0.0;
  const double da_t = -8.5;
  const double ltau = std::log(tau);
  const double tau2 = tau * tau;
  int c = 0;
  for (int s = 0; s < np; ++s) {
    const double h = (hi[s] - lo[s]) / 2;
    const double mid = (hi[s] + lo[s]) / 2;
    for (int k = 0; k < nn; ++k, ++c) {
      const double pc = mid + h * gk_nodes[k];
      const double zc = (pc - m) / tau;
      p[c] = pc;
      z[c] = zc;
      lw[c] = std::log(gk_wk[k] * h) + R::dnorm(zc, 0.0, 1.0, 1) - ltau;
      const double tfrac = edge[s] ? (pc - a) / denom : 0.0;
      double gm = zc / tau, gt = (zc * zc - 1) / tau;
      if (edge[s]) {
        gm += -da_m / denom;
        gt += -da_t / denom;
      }
      const double dm = edge[s] ? (1 - tfrac) * da_m : 0.0;
      const double dt = edge[s] ? (1 - tfrac) * da_t : 0.0;
      const double dlphi = -zc / tau;
      dpsi_m[c] = dm;
      dpsi_t[c] = dt;
      glw_m[c] = gm + dlphi * dm;
      glw_t[c] = gt + dlphi * dt;
      const double dzm = (dm - 1) / tau;
      const double dzt = (dt - zc) / tau;
      const double d2zmt = -(dm - 1) / tau2;
      const double d2ztt = -2 * (dt - zc) / tau2;
      double amm = -dzm * dzm;
      double amt = -dzm * dzt - zc * d2zmt;
      double att = -dzt * dzt - zc * d2ztt + 1 / tau2;
      if (edge[s]) {
        amm -= (da_m * da_m) / (denom * denom);
        amt -= (da_m * da_t) / (denom * denom);
        att -= (da_t * da_t) / (denom * denom);
      }
      alw_mm[c] = amm;
      alw_mt[c] = amt;
      alw_tt[c] = att;
    }
  }
  // the closed upper tail, a node at +Inf
  const double zn = (xn - m) / tau;
  const double inf = R_PosInf;
  p[c] = inf;
  z[c] = inf;
  lw[c] = R::pnorm(zn, 0.0, 1.0, 0, 1);
  glw_m[c] = mr / tau;
  glw_t[c] = zn * mr / tau;
  dpsi_m[c] = 0.0;
  dpsi_t[c] = 0.0;
  const double lamp = mr * (mr - zn);
  alw_mm[c] = -lamp / tau2;
  alw_mt[c] = -lamp * zn / tau2 - mr / tau2;
  alw_tt[c] = -lamp * zn * zn / tau2 - 2 * mr * zn / tau2;
  return List::create(_["p"] = p, _["lw"] = lw, _["z"] = z,
                      _["glw_m"] = glw_m, _["glw_t"] = glw_t,
                      _["dpsi_m"] = dpsi_m, _["dpsi_t"] = dpsi_t,
                      _["alw_mm"] = alw_mm, _["alw_mt"] = alw_mt,
                      _["alw_tt"] = alw_tt);
}

// [[Rcpp::export]]
List marg_seg_shift_cpp(const NumericVector& xg, const NumericVector& p,
                        bool linear, double beta, double gamma, bool jseg,
                        double delta) {
  const int ng = xg.size(), C = p.size();
  NumericMatrix shift(ng, C), hinge(ng, C), step(ng, C), dsdpsi(ng, C);
  for (int c = 0; c < C; ++c) {
    const bool fin = R_finite(p[c]);
    for (int t = 0; t < ng; ++t) {
      const double x = xg[t];
      const double hg = fin ? std::max(x - p[c], 0.0) : 0.0;
      const double st = (fin && x >= p[c]) ? 1.0 : 0.0;
      hinge(t, c) = hg;
      step(t, c) = st;
      double s = gamma * hg;
      if (linear) s = beta * x + s;
      if (jseg) s += delta * st;
      shift(t, c) = s;
      dsdpsi(t, c) = hg > 0 ? -gamma : 0.0;
    }
  }
  return List::create(_["shift"] = shift, _["hinge"] = hinge,
                      _["step"] = step, _["dshift_dpsi"] = dsdpsi);
}

static double lse(const std::vector<double>& a) {
  double mx = R_NegInf;
  for (double v : a) if (v > mx) mx = v;
  if (!R_finite(mx)) return R_NegInf;
  double s = 0;
  for (double v : a) s += std::exp(v - mx);
  return mx + std::log(s);
}

// The forward accumulation of one group: the log-likelihood of each
// observation given the ones before it, and its derivatives in the term's
// own coefficients (`D`, one ng by C matrix per coefficient) and in (m, tau)
// through the node weights and the node motion.
// [[Rcpp::export]]
List marg_seg_forward_cpp(const NumericVector& lw, const NumericMatrix& LD,
                          const NumericMatrix& SC, const List& D,
                          const NumericMatrix& dsdpsi,
                          const NumericVector& glw_m,
                          const NumericVector& glw_t,
                          const NumericVector& dpsi_m,
                          const NumericVector& dpsi_t) {
  const int ng = LD.nrow(), C = LD.ncol(), no = D.size();
  NumericVector ll(ng), jm(ng), jt(ng);
  NumericMatrix jo(ng, no);
  std::vector<NumericMatrix> Dm;
  for (int q = 0; q < no; ++q) Dm.push_back(as<NumericMatrix>(D[q]));
  std::vector<double> A(lw.begin(), lw.end()), A2(C), u(C), w(C), accp(C, 0.0),
    ap2(C);
  std::vector<std::vector<double>> acc(no, std::vector<double>(C, 0.0)),
    a2(no, std::vector<double>(C));
  double tot = lse(A);
  for (int c = 0; c < C; ++c) u[c] = std::exp(A[c] - tot);
  for (int t = 0; t < ng; ++t) {
    for (int c = 0; c < C; ++c) A2[c] = A[c] + LD(t, c);
    const double tot2 = lse(A2);
    for (int c = 0; c < C; ++c) w[c] = std::exp(A2[c] - tot2);
    ll[t] = tot2 - tot;
    for (int q = 0; q < no; ++q) {
      double s1 = 0, s0 = 0;
      for (int c = 0; c < C; ++c) {
        a2[q][c] = acc[q][c] + SC(t, c) * Dm[q](t, c);
        s1 += w[c] * a2[q][c];
        s0 += u[c] * acc[q][c];
      }
      jo(t, q) = s1 - s0;
      acc[q].swap(a2[q]);
    }
    double sm1 = 0, sm0 = 0, st1 = 0, st0 = 0;
    for (int c = 0; c < C; ++c) {
      ap2[c] = accp[c] + SC(t, c) * dsdpsi(t, c);
      sm1 += w[c] * (glw_m[c] + ap2[c] * dpsi_m[c]);
      sm0 += u[c] * (glw_m[c] + accp[c] * dpsi_m[c]);
      st1 += w[c] * (glw_t[c] + ap2[c] * dpsi_t[c]);
      st0 += u[c] * (glw_t[c] + accp[c] * dpsi_t[c]);
    }
    jm[t] = sm1 - sm0;
    jt[t] = st1 - st0;
    accp.swap(ap2);
    A.swap(A2);
    tot = tot2;
    u.swap(w);
  }
  return List::create(_["loglik"] = ll, _["own"] = jo, _["m"] = jm,
                      _["tau"] = jt);
}
