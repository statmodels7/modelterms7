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

// The general recursion of the submodel route: every coefficient of the
// filter read per observation, so a parameter carrying a subformula varies
// with t and the scalar kernel above does not apply. The derivative rows
// arrive already restricted to each group's active coordinates by the
// R-side preparation (one vectorized scatter per parameter per group), so
// the loop here allocates nothing per step; the callbacks into R remain,
// for the reason stated above the scalar kernel.
//
// Each element of `groups` carries: `rows` (1-based observation indices in
// time order), `act` (1-based positions of the group's active coordinates
// among the term's np), `dom` (m x na, the derivative row of the level at
// each observation), `da` (p matrices m x na, the same for each loading),
// `dbl` (q matrices m x na, the derivative row of each autoregressive
// coefficient, already chained through Levinson-Durbin), and `f0`/`df0`
// (the starting level and its derivative, from the group's first
// observation). `om`, `A` and `B` are the per-observation values of the
// level, the loadings and the autoregressive coefficients.
//
// [[Rcpp::export]]
List gas_filter_sub_cpp(NumericVector eta, List groups, int p, int q,
                        NumericVector om, NumericMatrix A, NumericMatrix B,
                        int np, Function score, Function curvature) {
    int n = eta.size();
    NumericVector eta_out(n);
    NumericMatrix jac(n, np);

    for (int g = 0; g < groups.size(); ++g) {
        List grp = groups[g];
        IntegerVector rows = grp["rows"];
        IntegerVector act = grp["act"];
        NumericMatrix dom = grp["dom"];
        List da = grp["da"];
        List dbl = grp["dbl"];
        double f0 = as<double>(grp["f0"]);
        NumericVector df0 = grp["df0"];
        int m = rows.size();
        int na = act.size();

        std::vector<double> f(m, 0.0), s(m, 0.0);
        std::vector<double> df(m * na, 0.0), ds(m * na, 0.0);
        std::vector<double> dft(na, 0.0);

        for (int t = 0; t < m; ++t) {
            int r = rows[t];
            double ft = om[r - 1];
            for (int k = 0; k < na; ++k) dft[k] = dom(t, k);

            for (int i = 1; i <= p; ++i) {
                double a_r = A(r - 1, i - 1);
                NumericMatrix dai = da[i - 1];
                double s_lag = (t - i >= 0) ? s[t - i] : 0.0;
                ft += a_r * s_lag;
                if (t - i >= 0) {
                    for (int k = 0; k < na; ++k) {
                        dft[k] += a_r * ds[(t - i) * na + k];
                    }
                }
                for (int k = 0; k < na; ++k) {
                    dft[k] += s_lag * dai(t, k);
                }
            }

            for (int j = 1; j <= q; ++j) {
                double b_r = B(r - 1, j - 1);
                NumericMatrix dbj = dbl[j - 1];
                double f_lag = (t - j >= 0) ? f[t - j] : f0;
                ft += b_r * f_lag;
                for (int k = 0; k < na; ++k) {
                    double df_lag = (t - j >= 0) ? df[(t - j) * na + k]
                                                 : df0[k];
                    dft[k] += b_r * df_lag + f_lag * dbj(t, k);
                }
            }

            f[t] = ft;
            for (int k = 0; k < na; ++k) df[t * na + k] = dft[k];

            double e_t = eta[r - 1] + ft;
            double s_t = as<double>(score(e_t, r));
            double c_t = as<double>(curvature(e_t, r));
            s[t] = s_t;
            for (int k = 0; k < na; ++k) ds[t * na + k] = c_t * dft[k];

            eta_out[r - 1] = e_t;
            for (int k = 0; k < na; ++k) jac(r - 1, act[k] - 1) = dft[k];
        }
    }

    return List::create(_["eta"] = eta_out, _["jacobian"] = jac);
}
