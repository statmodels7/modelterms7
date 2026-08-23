# The Derivative of a Filtered Predictor in the Static One

How the predictor a structural term produces moves when the static part
of the predictor moves, one row per observation and one column per
direction the caller supplies.

## Usage

``` r
term_static_deriv(term, curv, X, psi, ...)
```

## Arguments

- term:

  A built structural term.

- curv:

  The curvature at each predictor, as `term_filter` returns it.

- X:

  The directions to propagate, one column each – ordinarily the
  equation's design.

- psi:

  The term's parameters, on the parameter scale.

- ...:

  Passed to methods.

## Value

A matrix of `X`'s dimensions, or `NULL`.

## Details

A score-driven term's level is driven by scores read AT the predictor
the recursion is producing, so a coefficient in the same equation
reaches the level as well as the static part: writing \\f_t\\ for the
level and \\x_t\\ for a row of the design, \$\$\frac{\partial
\eta_t}{\partial \beta} = x_t + \frac{\partial f_t}{\partial \beta},\$\$
and the second piece obeys the recursion the filter already runs,
\$\$\frac{\partial f_t}{\partial \beta} = \sum_i \alpha_i \\
\ell''\_{t-i} \frac{\partial \eta\_{t-i}}{\partial \beta} + \sum_j
\beta_j \frac{\partial f\_{t-j}}{\partial \beta}.\$\$ The curvature it
needs is the one
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
returns, so no callback is evaluated here and the pass is arithmetic
alone.

Without it a standard error of the predictor counts the static part
only. Measured on a score-driven mean with one covariate beside it, that
understates the standard error by about a quarter.

The base method returns `NULL`: a term that is not a filter carries no
state, so the derivative is the design row itself and the caller needs
nothing from the term.

## See also

[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md),
[`term_adjoint`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20), x = rnorm(20))
term <- term_build(gas(p = 1, q = 1, time = t), dd)
psi <- list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5)
out <- term_filter(term, eta = rep(0, 20), y = dd$y,
                   score = function(e, i) dd$y[i] - e,
                   curvature = function(e, i) -1, psi = psi)
D <- term_static_deriv(term, out$curv, cbind(1, dd$x), psi)
dim(D)
#> [1] 20  2
```
