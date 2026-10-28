# The Filter's Parameters at Rows Outside the Fitting Data

Evaluates the recursion's coefficients \\(\omega_t, \kappa_t, \xi_t)\\
at each row of `newdata`, so that
[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)
can carry a score-driven filter past the end of its series. A parameter
carrying a subformula is read by calling
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
on each of its sub-terms, which reapplies the `build_info` that each
sub-term recorded at build time; a basis or a set of contrasts is
therefore not relearned on the new rows. A parameter with no subformula
is constant over the rows.

## Usage

``` r
.gas_sub_new(term, u, newdata)
```

## Arguments

- term:

  A built `GasTerm` whose `build_info` carries at least one subformula.
  A term without a subformula is handled by `.gas_coefs()`.

- u:

  The term's parameter vector, in
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  order, on the mixed scale described above.

- newdata:

  A data frame of the rows at which to read the parameters, carrying
  every variable that the subformulas name.

## Value

A list with the following elements:

- `omega`:

  numeric, one \\\omega_t\\ per row of `newdata`.

- `A`:

  `nrow(newdata)` by `p`, the score loadings.

- `B`:

  `nrow(newdata)` by `q`, the autoregressive coefficients, obtained by
  applying Levinson-Durbin row by row. At \\q = 0\\ it is a zero-column
  matrix.

- `d_omega`, `d_A`, `d_B`:

  the derivatives of `omega`, of each column of `A` and of each column
  of `B` with respect to `u`, each a `nrow(newdata)` by `length(u)`
  matrix (`d_A` and `d_B` are lists of them).

## Details

For a developed parameter \\j\\ the value at row \\t\\ is \\\psi\_{j,t}
= g_j^{-1}(z_t' \gamma_j)\\, where \\g_j\\ is the parameter's own link
and \\\gamma_j\\ its coefficients. The partial autocorrelations obtained
in this way are mapped to autoregressive coefficients by the
Levinson-Durbin recursion, so `B` holds autoregressive coefficients
where `u` holds the coordinates of the partial autocorrelations.

The vector `u` has the two scales of the parameter vector that
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
receives. An entry for a parameter with no subformula is the parameter
itself, on the parameter scale. The entries of a developed parameter are
its coefficients on the unconstrained scale, and the chart is applied to
them here. For example, with \\p = q = 1\\, the development
`omega ~ 0 + z` and scalar \\\kappa_1\\ and partial autocorrelation,
`u = c(0.3, 0.4, 0.7)` gives \\\omega_t = 0.3 z_t\\, \\\kappa_1 = 0.4\\
and a partial autocorrelation of 0.7.

## See also

[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md),
the only caller;
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
for the same quantities inside the series.
