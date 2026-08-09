#include <Rcpp.h>
using namespace Rcpp;

// The normalized forward recursion of a Markov regime term, and the
// derivative propagated beside it.
//
// Unlike the score-driven filter of src/gas_filter.cpp, nothing here has
// to call back into R: a regime shifts the predictor by a level of its
// own, so the log-density and the score of every observation under every
// regime are known before the recursion starts. They arrive as the two
// n by k matrices LF and SC, computed by k vectorized calls.
//
// Written out, with a the filtered state and da its derivative:
//   w_j    = exp(lf_j - max lf)
//   pred   = a P                       (the stationary vector at t = 1)
//   atil_j = w_j pred_j,   c = sum atil
//   loglik = log c + max lf,   jac = (d c) / c
//   a      = atil / c,   da = (d atil - (d c) a') / c
//
// [[Rcpp::export]]
List regime_forward_cpp(const List& order, const NumericMatrix& LF,
                        const NumericMatrix& SC, const NumericMatrix& dmu,
                        const NumericMatrix& P, const List& dP,
                        const NumericVector& delta,
                        const NumericMatrix& ddelta) {
  const int n = LF.nrow(), k = LF.ncol(), np = dmu.ncol();
  NumericVector loglik(n);
  NumericMatrix jac(n, np);
  if (n == 0) return List::create(_["loglik"] = loglik, _["jacobian"] = jac);

  // the chain's derivatives flattened once, so the inner loop indexes
  // contiguous memory rather than an R list
  std::vector<double> dPf((size_t)np * k * k, 0.0);
  for (int i = 0; i < np; i++) {
    NumericMatrix Di = dP[i];
    for (int c = 0; c < k; c++)
      for (int r = 0; r < k; r++)
        dPf[((size_t)i * k + c) * k + r] = Di(r, c);
  }

  std::vector<double> a(k), pred(k), atil(k), w(k), sc(k);
  std::vector<double> da((size_t)np * k), dpred((size_t)np * k),
                      datil((size_t)np * k), dw((size_t)np * k), dct(np);

  for (int gi = 0; gi < order.size(); gi++) {
    IntegerVector rows = order[gi];
    const int m = rows.size();
    for (int j = 0; j < k; j++) a[j] = delta[j];
    for (int i = 0; i < np; i++)
      for (int j = 0; j < k; j++) da[(size_t)i * k + j] = ddelta(i, j);

    for (int t = 0; t < m; t++) {
      const int row = rows[t] - 1;

      double mx = LF(row, 0);
      for (int j = 1; j < k; j++) if (LF(row, j) > mx) mx = LF(row, j);
      for (int j = 0; j < k; j++) {
        w[j] = std::exp(LF(row, j) - mx);
        sc[j] = SC(row, j);
      }
      // dw(i, j) = dmu(j, i) * w_j * s_j
      for (int i = 0; i < np; i++)
        for (int j = 0; j < k; j++)
          dw[(size_t)i * k + j] = dmu(j, i) * w[j] * sc[j];

      if (t == 0) {
        for (int j = 0; j < k; j++) pred[j] = a[j];
        for (size_t z = 0; z < (size_t)np * k; z++) dpred[z] = da[z];
      } else {
        for (int j = 0; j < k; j++) {
          double s = 0.0;
          for (int r = 0; r < k; r++) s += a[r] * P(r, j);
          pred[j] = s;
        }
        for (int i = 0; i < np; i++) {
          const double* dai = &da[(size_t)i * k];
          const double* Di = &dPf[(size_t)i * k * k];
          for (int j = 0; j < k; j++) {
            double s = 0.0;
            for (int r = 0; r < k; r++)
              s += dai[r] * P(r, j) + a[r] * Di[(size_t)j * k + r];
            dpred[(size_t)i * k + j] = s;
          }
        }
      }

      double ct = 0.0;
      for (int j = 0; j < k; j++) { atil[j] = w[j] * pred[j]; ct += atil[j]; }
      for (int i = 0; i < np; i++) {
        double s = 0.0;
        for (int j = 0; j < k; j++) {
          const size_t z = (size_t)i * k + j;
          datil[z] = dw[z] * pred[j] + dpred[z] * w[j];
          s += datil[z];
        }
        dct[i] = s;
      }

      loglik[row] = std::log(ct) + mx;
      for (int i = 0; i < np; i++) jac(row, i) = dct[i] / ct;

      for (int j = 0; j < k; j++) a[j] = atil[j] / ct;
      for (int i = 0; i < np; i++)
        for (int j = 0; j < k; j++)
          da[(size_t)i * k + j] = (datil[(size_t)i * k + j] - dct[i] * a[j]) / ct;
    }
  }

  return List::create(_["loglik"] = loglik, _["jacobian"] = jac);
}
