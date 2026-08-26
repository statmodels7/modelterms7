# The Filter's Parameters at Rows Outside the Fitting Data

Evaluates the recursion's coefficients \\(\omega_t, \alpha_t, \beta_t)\\
at each row of `newdata`, so that
[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)
can carry a score-driven filter past the end of its series. A parameter
carrying a subformula is read by calling
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
on each of its sub-terms, which reapplies the blueprint each one
recorded at build time; a basis or a set of contrasts is therefore not
relearned on the new rows. A parameter with no subformula is constant
along them.

## Usage

``` r
.gas_sub_new(term, u, newdata)
```

## Arguments

- term:

  A built `GasTerm` whose blueprint carries at least one subformula.
  With no development the caller uses `.gas_coefs()` instead.

- u:

  The term's parameter vector, in
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  order, on the mixed scale above.

- newdata:

  A data frame of the rows to read at, carrying every variable the
  subformulas name.

## Value

A list of three:

- `om`:

  numeric, one \\\omega_t\\ per row of `newdata`.

- `A`:

  `nrow(newdata)` by `p`, the score loadings.

- `B`:

  `nrow(newdata)` by `q`, the autoregressive coefficients,
  Levinson-Durbin applied row by row. A zero-column matrix at \\q = 0\\.

## Details

For a developed parameter \\j\\ the value at row \\t\\ is \\\psi\_{j,t}
= g_j^{-1}(z_t' \gamma_j)\\ with \\g_j\\ the parameter's own link and
\\\gamma_j\\ its coefficients, and the persistence goes on to
Levinson-Durbin, so `B` holds autoregressive coefficients where `u`
holds partial autocorrelations.

`u` mixes two scales, matching what
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
is handed: a parameter with no subformula is taken as it stands, on the
parameter scale, while a developed parameter's entries are its
coefficients on the unconstrained scale and go through the chart here.
With \\p = q = 1\\ and no development, `u = c(0.3, 0.4, 0.7)` gives
\\\omega = 0.3\\, \\\alpha_1 = 0.4\\ and a partial autocorrelation of
0.7.

## See also

[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md),
the only caller;
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
for the same quantities inside the series.
