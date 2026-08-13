# Second Derivatives of a Structural Term's Predictor

The Jacobian of the predictor the term produces with respect to a
caller's unknowns, and the second derivative of that predictor
contracted against a caller's weights. It is what an observed
information needs and
[`term_adjoint`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
does not give.

## Usage

``` r
term_curvature(term, eta, y, score, curvature, psi, g, seed, blocks, ...)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static part of the predictor.

- y:

  The response.

- score, curvature:

  The callbacks of
  [`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The term's parameters, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The weights the second derivative is contracted against, one per
  observation.

- seed:

  The derivative of the static predictor in the caller's unknowns, one
  row per observation.

- blocks:

  A function of the predictor, the index and the current Jacobian row,
  returning `cross` and `M`.

- ...:

  Passed to methods.

## Value

A list with `jacobian`, the derivative of the predictor in the caller's
unknowns, and `curvature`, the contracted second derivative.

## Details

The gradient of a model carrying a filter needs no second derivative and
no forward Jacobian: the reverse recursion answers it in one pass. The
curvature needs both. Writing \\u\\ for the caller's unknowns, \\D_t =
\partial e_t/\partial u\\ and \\E_t = \partial^2 e_t/\partial u \partial
u^\top\\, the observed information of the model is

\$\$\sum_t w_t \sum\_{q,r} \ell\_{qr,t} V\_{q,t}^\top V\_{r,t} + \sum_t
w_t \ell\_{p,t} E_t,\$\$

whose first sum is a weighted crossproduct the caller assembles and
whose second is what this returns, the weights \\g_t = w_t\ell\_{p,t}\\
being supplied.

`seed` is \\\partial \eta^{0}/\partial u\\, one row per observation: the
caller says how its unknowns reach the static predictor and the term
knows nothing else about them. `blocks` is where the model's own
derivatives enter, since the score the recursion is driven by depends on
every equation and not only on the predictor it is read at. It is called
once per observation with the predictor there and the Jacobian the
recursion has reached, and returns the two quantities that seed the
first and second derivatives of that score,

\$\$\dot S_t = \ell\_{pp,t}D_t + \texttt{cross}, \qquad \ddot S_t =
\ell\_{pp,t}E_t + \texttt{M},\$\$

with `cross` \\= \sum\_{q\ne p}\ell\_{pq,t}C\_{q,t}\\ and `M` \\=
\sum\_{r,r'}\ell\_{prr',t}V\_{r,t}^\top V\_{r',t}\\, the third
derivatives of the log-density in the predictors. A model of one
equation has `cross` zero and `M` equal to \\\ell\_{ppp}D^\top D\\.

## See also

[`term_adjoint`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md),
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20))
term <- term_build(gas(p = 1, q = 1, time = t), dd)

# a gaussian mean of unit variance and one equation, so the model
# contributes no cross term and no third derivative
m <- 3L
out <- term_curvature(
  term, eta = rep(0, 20), y = dd$y,
  score = function(e, i) dd$y[i] - e,
  curvature = function(e, i) -1,
  psi = list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5),
  g = rep(1, 20), seed = matrix(0, 20, m),
  blocks = function(e, i, D) list(cross = numeric(m),
                                  M = matrix(0, m, m)))
dim(out$jacobian)
#> [1] 20  3
dim(out$curvature)
#> [1] 3 3
```
