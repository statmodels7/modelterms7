# Second Derivatives of a Structural Term's Predictor

The Jacobian of the predictor the term produces with respect to a
caller's unknowns, and the second derivative of that predictor
contracted against a caller's weights. This second derivative is the
piece of the Hessian of the log-likelihood that
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
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
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The term's parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The weights the second derivative is contracted against, one per
  observation.

- seed:

  The derivative of the static predictor in the caller's unknowns, one
  row per observation.

- blocks:

  A function `function(e, i, D, act)` of the predictor, the observation
  index, the current Jacobian row restricted to the active set and that
  set, returning `cross` and `M` on the active set. A function of the
  first three arguments only is accepted and is given the full row.

- ...:

  Passed to methods.

## Value

A list of two: `jacobian`, an `n` by `ncol(seed)` matrix holding the
derivative of the predictor in the caller's unknowns, and `curvature`,
the `ncol(seed)` by `ncol(seed)` second derivative contracted against
`g`. Both are exact.

## Details

The gradient of a model carrying a filter needs no second derivative and
no forward Jacobian: the reverse recursion gives it in one pass. The
curvature needs both. Writing \\u\\ for the caller's unknowns, \\p\\ for
the distribution parameter whose predictor carries the filter, \\D_t =
\partial e_t/\partial u\\ for the derivative of that predictor,
\\V\_{q,t}\\ for the derivative of the predictor of parameter \\q\\
(with \\V\_{p,t} = D_t\\) and \\E_t = \partial^2 e_t/\partial u \partial
u^\top\\, the Hessian of the log-likelihood in \\u\\ is

\$\$\sum_t w_t \sum\_{q,r} \ell\_{qr,t} V\_{q,t}^\top V\_{r,t} + \sum_t
w_t \ell\_{p,t} E_t,\$\$

and the observed information is its negative. The first sum is a
weighted crossproduct that the caller assembles, and the second is what
this returns, the weights \\g_t = w_t\ell\_{p,t}\\ being supplied.

`seed` is \\\partial \eta^{0}/\partial u\\, one row per observation: the
caller supplies how its unknowns reach the static predictor, and the
term holds no other information about them. `blocks` is where the
model's own derivatives enter, since the score that drives the recursion
depends on every equation, beyond the predictor at which it is read. It
is called once per observation, with the predictor there, the index, the
Jacobian row that the recursion has reached (restricted to the active
coordinates of the group) and the indices of those coordinates, and it
returns the two quantities that seed the first and second derivatives of
that score,

\$\$\dot S_t = \ell\_{pp,t}D_t + \texttt{cross}, \qquad \ddot S_t =
\ell\_{pp,t}E_t + \texttt{M},\$\$

with `cross` \\= \sum\_{q\ne p}\ell\_{pq,t}V\_{q,t}\\ and `M` \\=
\sum\_{r,r'}\ell\_{prr',t}V\_{r,t}^\top V\_{r',t}\\, the third
derivatives of the log-density in the predictors. A model of one
equation has `cross` zero and `M` equal to \\\ell\_{ppp}D^\top D\\.

## See also

[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
for the first order,
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
for the third,
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
for the callbacks.

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
  psi = list(omega = 0.1, kappa1 = 0.3, pacf1 = 0.5),
  g = rep(1, 20), seed = matrix(0, 20, m),
  blocks = function(e, i, D) list(cross = numeric(m),
                                  M = matrix(0, m, m)))
dim(out$jacobian)
#> [1] 20  3
dim(out$curvature)
#> [1] 3 3
```
