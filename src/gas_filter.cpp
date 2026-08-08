#include <Rcpp.h>
using namespace Rcpp;

// The score-driven recursion and the derivative of its state, propagated
// together. The two callbacks into R remain -- the score and the curvature
// belong to whatever distribution the model carries, and are functions the
// caller supplies -- but they were measured at 17 to 27 per cent of the R
// loop's time, so the arithmetic around them is what this removes.
//
// Every index arriving from R is 1-based and is used to address the caller's
// rows, which is why `row` is passed to the callbacks unchanged.
//
// [[Rcpp::export]]
List gas_filter_cpp(NumericVector eta, List order, int p, int q,
                    double omega, NumericVector a, NumericVector b,
                    NumericMatrix db, double f0, NumericVector df0,
                    IntegerVector i_a, int np,
                    Function score, Function curvature) {
    int n = eta.size();
    NumericVector eta_out(n);
    NumericMatrix jac(n, np);

    for (int g = 0; g < order.size(); ++g) {
        IntegerVector rows = order[g];
        int m = rows.size();
        std::vector<double> f(m, 0.0), s(m, 0.0);
        std::vector<double> df(m * np, 0.0), ds(m * np, 0.0);
        std::vector<double> dft(np, 0.0);

        for (int t = 0; t < m; ++t) {
            double ft = omega;
            std::fill(dft.begin(), dft.end(), 0.0);
            dft[0] = 1.0;

            for (int i = 1; i <= p; ++i) {
                double s_lag = (t - i >= 0) ? s[t - i] : 0.0;
                ft += a[i - 1] * s_lag;
                if (t - i >= 0) {
                    for (int k = 0; k < np; ++k) {
                        dft[k] += a[i - 1] * ds[(t - i) * np + k];
                    }
                }
                dft[i_a[i - 1] - 1] += s_lag;
            }

            for (int j = 1; j <= q; ++j) {
                double f_lag = (t - j >= 0) ? f[t - j] : f0;
                ft += b[j - 1] * f_lag;
                for (int k = 0; k < np; ++k) {
                    double df_lag = (t - j >= 0) ? df[(t - j) * np + k] : df0[k];
                    dft[k] += b[j - 1] * df_lag + db(j - 1, k) * f_lag;
                }
            }

            f[t] = ft;
            for (int k = 0; k < np; ++k) df[t * np + k] = dft[k];

            int row = rows[t];
            double e_t = eta[row - 1] + ft;
            double s_t = as<double>(score(e_t, row));
            double c_t = as<double>(curvature(e_t, row));
            s[t] = s_t;
            for (int k = 0; k < np; ++k) ds[t * np + k] = c_t * dft[k];

            eta_out[row - 1] = e_t;
            for (int k = 0; k < np; ++k) jac(row - 1, k) = dft[k];
        }
    }

    return List::create(_["eta"] = eta_out, _["jacobian"] = jac);
}
