#include <Rcpp.h>
using namespace Rcpp;

// The working block of a segmented, stepmented or joint term, and the
// contribution it linearizes. Written as one pass per column into the
// result matrix: the R form of the same arithmetic allocates six vectors
// per break-point, which at n = 1e5 dominates the elementwise work.
//
// kind: 0 = seg (continuous), 1 = jump (discontinuous), 2 = both.
// Column order matches .seg_names(): the linear column if present, then
// the slope changes (seg, jseg), then either the break-point Jacobians
// (seg) or the kappa columns followed by the g columns (jump, jseg).
//
// [[Rcpp::export]]
List seg_block_cpp(int kind, NumericVector xv, NumericVector psi,
                   NumericVector del, NumericVector kap, double lin,
                   bool linear, double floor_w) {
  const R_xlen_t n = xv.size();
  const int npsi = psi.size();
  const bool has_delta = (kind == 0 || kind == 2);
  const bool has_jump = (kind == 1 || kind == 2);

  int ncol = (linear ? 1 : 0);
  if (has_delta) ncol += npsi;
  ncol += (kind == 0) ? npsi : 2 * npsi;

  NumericMatrix X(n, ncol);
  NumericVector value(n);
  if (n == 0) return List::create(_["X"] = X, _["value"] = value);

  double* v = REAL(value);
  const double* x = REAL(xv);
  double* Xp = REAL(X);

  if (linear) {
    double* c0 = Xp;
    for (R_xlen_t i = 0; i < n; i++) {
      c0[i] = x[i];
      v[i] = lin * x[i];
    }
  }
  int col = linear ? 1 : 0;

  if (has_delta) {
    for (int j = 0; j < npsi; j++) {
      const double p = psi[j], d = del[j];
      double* cj = Xp + (R_xlen_t)(col + j) * n;
      for (R_xlen_t i = 0; i < n; i++) {
        const double t = x[i] - p;
        const double u = (t > 0.0) ? t : 0.0;
        cj[i] = u;
        v[i] += d * u;
      }
    }
    col += npsi;
  }

  if (kind == 0) {
    for (int j = 0; j < npsi; j++) {
      const double p = psi[j], d = del[j];
      double* cj = Xp + (R_xlen_t)(col + j) * n;
      for (R_xlen_t i = 0; i < n; i++) cj[i] = (x[i] > p) ? -d : 0.0;
    }
  } else if (has_jump) {
    for (int j = 0; j < npsi; j++) {
      const double p = psi[j], k = kap[j];
      double* cz = Xp + (R_xlen_t)(col + j) * n;
      double* cw = Xp + (R_xlen_t)(col + npsi + j) * n;
      for (R_xlen_t i = 0; i < n; i++) {
        double den = 2.0 * std::fabs(x[i] - p);
        if (den < floor_w) den = floor_w;
        const double w = 1.0 / den;
        cw[i] = w;
        cz[i] = x[i] * w + 0.5;
        if (x[i] > p) v[i] += k;
      }
    }
  }

  return List::create(_["X"] = X, _["value"] = value);
}
