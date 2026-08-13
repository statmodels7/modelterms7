# Differentiate a Structural Term Backwards

The derivative of a caller's objective with respect to the static
predictor the term was handed, and with respect to the sequence of
scores it was given, both accounting for the recursion.

## Usage

``` r
term_adjoint(term, eta, y, score, curvature, psi, g, ...)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static part of the predictor, as
  [`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
  takes it.

- y:

  The response.

- score, curvature:

  The callbacks of
  [`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The term's parameters, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The direct derivative of the objective in the predictor the term
  produced, one value per observation.

- ...:

  Passed to methods.

## Value

A list with `deta` and `dscore`, each one value per observation.

## Details

[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
returns the derivative of the predictor in the term's OWN parameters,
which is what estimating those needs. It is not what estimating the
coefficients of the same equation needs: the level at one time is driven
by the scores at earlier ones, which are read at predictors those
coefficients also enter, so the derivative of the predictor in a
coefficient carries a term the block does not. Propagating that forward
would cost one derivative array per coefficient; the reverse recursion
here costs one pass whatever their number, and is exact.

Two derivatives are returned rather than one because the score the
caller supplies depends on more than the predictor it is read at.
Writing \\s_t\\ for that score and \\\bar{s}\_t\\ for `dscore`, the
derivative of the objective in anything the score depends on is

\$\$\frac{\partial L}{\partial \theta} = \left.\frac{\partial
L}{\partial \theta}\right\|\_{\mathrm{direct}} + \sum_t \bar{s}\_t
\frac{\partial s_t}{\partial \theta},\$\$

so a model layer whose score is the derivative of its log-likelihood in
one distribution parameter obtains the derivative in the predictor of
ANOTHER by multiplying `dscore` by the mixed second derivative of that
log-likelihood. `deta` is that formula applied to the term's own
equation, where the second factor is the curvature.

## See also

[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20))
term <- term_build(gas(p = 1, q = 1, time = t), dd)
out <- term_adjoint(term, eta = rep(0, 20), y = dd$y,
                    score = function(e, i) dd$y[i] - e,
                    curvature = function(e, i) -1,
                    psi = list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5),
                    g = rep(1, 20))
head(out$deta, 3)
#> [1] 0.625 0.625 0.625
```
