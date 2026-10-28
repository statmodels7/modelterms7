# The Derivative of a Filtered Predictor in the Static One

Returns how the predictor that a structural term produces moves when the
static part of the predictor moves, one row per observation and one
column per direction that the caller supplies.

## Usage

``` r
term_static_deriv(term, curvature, X, psi, ...)
```

## Arguments

- term:

  A built structural term.

- curvature:

  The curvature at each predictor, as `term_filter` returns it.

- X:

  The directions to propagate, one column each. Ordinarily `X` is the
  equation's design, so the result has one row per observation and one
  column per coefficient.

- psi:

  The term's parameters, on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A matrix of `X`'s dimensions, holding the full derivative \\\partial
\eta_t/\partial \beta = x_t + \partial f_t/\partial \beta\\ for each
column of `X`, so that the design is already included. The value is
`NULL` for a term that carries no state.

## Details

A score-driven term's level is driven by scores read at the predictor
that the recursion is producing, so a coefficient in the same equation
reaches the level as well as the static part: writing \\f_t\\ for the
level and \\x_t\\ for a row of the design, \$\$\frac{\partial
\eta_t}{\partial \beta} = x_t + \frac{\partial f_t}{\partial \beta},\$\$
and the second piece obeys the recursion that the filter already runs,
\$\$\frac{\partial f_t}{\partial \beta} = \sum_i \kappa_i \\
\ell''\_{t-i} \frac{\partial \eta\_{t-i}}{\partial \beta} + \sum_j \xi_j
\frac{\partial f\_{t-j}}{\partial \beta},\$\$ with \\\kappa_i\\ the
score loadings, \\\xi_j\\ the autoregressive coefficients and
\\\ell''\_t\\ the curvature. The curvature is the one that
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
returns, so this pass does not evaluate any callback and consists only
of arithmetic.

Without this derivative, the standard error of the predictor accounts
for the static part only and leaves out the part that runs through the
filter.

The base method returns `NULL`. A term that is not a filter carries no
state, so the derivative is the design row itself and the caller needs
nothing from the term.

## See also

[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md),
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20), x = rnorm(20))
term <- term_build(gas(p = 1, q = 1, time = t), dd)
psi <- list(omega = 0.1, kappa1 = 0.3, pacf1 = 0.5)
out <- term_filter(term, eta = rep(0, 20), y = dd$y,
                   score = function(e, i) dd$y[i] - e,
                   curvature = function(e, i) -1, psi = psi)
D <- term_static_deriv(term, out$curv, cbind(1, dd$x), psi)
dim(D)
#> [1] 20  2
```
