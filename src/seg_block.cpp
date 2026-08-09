#include <Rcpp.h>
using namespace Rcpp;

// The working block of a segmented, stepmented or joint term, and the
// contribution it linearizes. Written as one pass per column into the
// result matrix: the R form of the same arithmetic allocates several
// vectors per break-point, which at n = 1e5 dominates the elementwise
// work.
//
// kind: 0 = seg (continuous), 1 = jump (discontinuous), 2 = both.
// Column order matches .seg_names(): the linear column if present, then
// the slope changes (seg, jseg), then either the break-point Jacobians
// (seg) or the kappa columns followed by the g columns (jump, jseg).
//
// For a discontinuous term the indicator columns are built on the
// covariate RESCALED away from the break-point, after Fasola, Muggeo and
// Kuchenhoff (2018): with a scaling factor c the two intervals
// [lo, psi] and (psi, hi] are mapped onto [lo, psi - c(psi - lo)] and
// (psi + c(hi - psi), hi], which leaves a gap around psi and bounds the
// weight 1/(2|x' - psi|) without altering the model. The truncated line
// and the linear column stay on the original covariate, as the
// additional covariates of the paper's equation (14) do.
//
// [[Rcpp::export]]
List seg_block_cpp(int kind, NumericVector xv, NumericVector psi,
                   NumericVector del, NumericVector kap, NumericVector cvec,
                   double lin, bool linear, double lo, double hi) {
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

  // The truncated line and the linear column stay on the ORIGINAL
  // covariate, as the additional covariates of Fasola et al.'s equation
  // (14) do: computing them on the rescaled one was measured and is
  // worse, the displacement biasing the slope it is not needed for.
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
      const double p = psi[j], k = kap[j], c = cvec[j];
      const double keep = 1.0 - c;
      const double pp = p + c * (hi - p);
      double* cz = Xp + (R_xlen_t)(col + j) * n;
      double* cw = Xp + (R_xlen_t)(col + npsi + j) * n;
      for (R_xlen_t i = 0; i < n; i++) {
        const double xs = (x[i] <= p) ? lo + (x[i] - lo) * keep
                                      : pp + (x[i] - p) * keep;
        const double w = 1.0 / (2.0 * std::fabs(xs - p));
        cw[i] = w;
        cz[i] = xs * w + 0.5;
        if (x[i] > p) v[i] += k;
      }
    }
  }

  return List::create(_["X"] = X, _["value"] = value);
}
