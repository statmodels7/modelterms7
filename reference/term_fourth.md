# Fourth Derivatives of a Structural Term's Predictor, in Two Directions

[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
differentiated once more along a second direction. It is what the exact
HESSIAN of a marginal criterion needs when the model carries a term that
bends the predictor.

## Usage

``` r
term_fourth(
  term,
  eta,
  y,
  score,
  curvature,
  psi,
  g,
  seed,
  blocks,
  directions,
  ...
)
```

## Arguments

- term:

  A built term.

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

  The weights the fourth derivative is contracted against, one per
  observation.

- seed:

  The derivative of the static predictor in the caller's unknowns.

- blocks:

  A function of the predictor, the index, the current Jacobian row and
  the active set, returning `cross`, `M`, `dcurv`, `N`, `Q`, `P` and
  `cppp`.

- directions:

  A list of two numeric vectors, the directions \\v\\ and \\w\\ in the
  caller's unknowns. The result is symmetric in the two.

- ...:

  Passed to methods.

## Value

A list with `jacobian`, the derivative of the predictor in the caller's
unknowns; `dphi`, a list of two, the second derivative contracted
against each direction, one row per observation; `dpsi`, the third
derivative contracted against both, one row per observation; and
`curvature`, the fourth derivative contracted against `g` and both
directions.

## Details

A marginal criterion carries \\-\frac{1}{2}\log\|K\|\\ at the penalized
mode, so its gradient in one hyperparameter reads \\\partial K/\partial
u\\ along the direction the mode moves in, which is
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md),
and its second derivative in a PAIR reads the same object along two
directions. Where the predictor is \\X\beta\\ that is the family's
derivatives against the design and nothing else; where a filter produces
it, \\K\\ carries \\\sum_t w_t \ell\_{p,t}E_t\\ as well, and
differentiating that twice asks for \\\partial^4 e_t/\partial u^4\\.

**Only contracted, and only twice.** The full fourth derivative is an
\\m^4\\ array per observation and is never formed: what is propagated is
a matrix per observation, the same size as the curvature and therefore
the same \\O(nm^2)\\. A caller wanting a Hessian over \\n_h\\
hyperparameters calls this once per pair.

**Five states, not three.** The recursion carries \\F_t\\, \\\Phi_t\\,
\\\Psi^{(v)}\_t\\, \\\Psi^{(w)}\_t\\ and \\\Omega_t\\: BOTH third
derivatives, because the product rule at the fourth order pairs each
direction's third derivative with the other direction, so a recursion
holding one contraction cannot reach this order.

**What the model supplies.** Each order of differentiation pulls in one
more order of the family, so beyond the third order's `dcurv` and `N`
the `blocks` callback returns, at an observation,

\$\$\texttt{Q} = \sum\_{r,r'}\ell\_{pprr'}V_r^\top V\_{r'}, \qquad
\texttt{P} = \sum\_{r,r'}\Big(\sum\_{s,s'}\ell\_{prr'ss'} (V_s\cdot
v)(V\_{s'}\cdot w)\Big)V_r^\top V\_{r'},\$\$

with `cppp` \\= \ell\_{ppp}\\ beside them, and `N` a list of two, one
per direction. `P` is the only place the family's FIFTH derivative
enters; `cppp` is a scalar no contraction recovers, the level's own
curvature moving by it.

The base method returns zeros, so an additive term is covered without
writing anything. A structural term refuses, for the reason
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
gives: zero there would be a false statement rather than a true one.

## See also

[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
for the third order this differentiates,
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
for the second,
[`statmodels7::reml()`](https://statmodels7.github.io/statmodels7/reference/reml.html)
for the criterion whose Hessian asks for it.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20))
term <- term_build(linpar(~t), dd)
# an additive term bends no predictor, so every order above the first is
# zero and the base method says so
out <- term_fourth(term, rep(0, 20), dd$y,
                   score = function(e, i) dd$y[i] - e,
                   curvature = function(e, i) -1, psi = list(),
                   g = rep(1, 20), seed = matrix(0, 20, 2),
                   blocks = function(e, i, D, act) NULL,
                   directions = list(c(1, 0), c(0, 1)))
all(out$curvature == 0)
#> [1] TRUE
```
