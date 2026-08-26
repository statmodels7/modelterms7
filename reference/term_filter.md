# Apply a Structural Term to a Linear Predictor

Runs the term's recursion over the data and returns the predictor it
produces, together with the derivative of that predictor with respect to
the term's own parameters. This is the operation that makes a structural
term structural: the predictor at one observation depends on the others,
so it cannot be written as a block of columns.

## Usage

``` r
term_filter(term, eta, y, score, curvature, psi, ...)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static part of the linear predictor, one value per observation.

- y:

  The response.

- score:

  A function of the predictor returning the derivative of the
  log-likelihood with respect to it, one value per observation.

- curvature:

  A function of the predictor returning the second derivative of the
  log-likelihood with respect to it.

- psi:

  The term's parameters, on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A list with `eta`, the predictor the term produces, `jacobian`, an `n`
by `length(psi)` matrix of its derivatives with respect to `psi`, and
`curv`, the value of `curvature` at each predictor. That last one is
returned because the recursion evaluates it anyway, so a consumer
running a second pass at the same point, as
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
does, reads it instead of evaluating the callback again.

## Details

Writing \\\eta_t^{0}\\ for the static predictor supplied in `eta` and
\\\psi\\ for the term's parameters, the filter returns the pair

\$\$\eta_t = \eta_t^{0} + f_t(\psi), \qquad J\_{tj} = \frac{\partial
\eta_t}{\partial \psi_j},\$\$

where \\f_t\\ is the term's own recursion, driven by `score` and
`curvature` evaluated at the predictor already produced. Both are read
at \\\eta_t\\, so `curvature` is the second derivative \\\partial^{2}
\ell_t / \partial \eta^{2}\\ and is negative at an ordinary observation.
**Its sign is load-bearing**: passing the information, which is its
negative, returns a predictor and a Jacobian that are internally
consistent and wrong, with no error to say so.

The derivative is returned because the recursion is the only place it
can be computed. Propagating it beside the state costs one extra vector
per parameter and is exact; a model layer differencing the filter
instead would pay one pass per parameter and inherit the error of the
difference.

## See also

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
and
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
for the two structural shapes,
[`term_static_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_static_deriv.md)
for the same recursion in the static predictor,
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
for the reverse pass,
[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)
for the recursion past the series.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20))
term <- term_build(gas(p = 1, q = 1, time = t), dd)

# the score and curvature a Gaussian mean would supply
out <- term_filter(term, eta = rep(0, 20), y = dd$y,
                   score = function(e, i) dd$y[i] - e,
                   curvature = function(e, i) -1,
                   psi = list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5))
head(out$eta, 3)
#> [1]  0.20000000 -0.04793614  0.14550577
dim(out$jacobian)
#> [1] 20  3
```
