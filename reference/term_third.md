# Third Derivatives of a Structural Term's Predictor, in One Direction

Returns the second derivative of
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
differentiated once more along a single direction, and the derivative of
the term's Jacobian along that same direction. The exact gradient of a
marginal criterion needs this when a penalty covers the term's own
parameters.

## Usage

``` r
term_third(
  term,
  eta,
  y,
  score,
  curvature,
  psi,
  g,
  seed,
  blocks,
  direction,
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

  The weights against which the third derivative is contracted, one per
  observation.

- seed:

  The derivative of the static predictor in the caller's unknowns.

- blocks:

  A function of the predictor, the index, the current Jacobian row and
  the active set, returning `cross`, `M`, `dcurv` and `N`.

- direction:

  The direction \\v\\, in the caller's unknowns.

- ...:

  Passed to methods.

## Value

A list with `jacobian`, the derivative of the predictor in the caller's
unknowns; `dphi`, the second derivative contracted against `direction`,
one row per observation; and `curvature`, the third derivative
contracted against both `g` and `direction`.

## Details

A marginal criterion carries \\-\frac{1}{2}\log\|K\|\\ at the penalized
mode, so its gradient in a hyperparameter needs
\\\mathrm{tr}(M\\\partial K/\partial u\\\[v\])\\, with \\v\\ the
direction in which the mode moves. Where the predictor is \\X\beta\\
that contraction is the family's third derivative against the design and
nothing else. Where a filter produces the predictor, \\K\\ carries
\\\sum_t w_t \ell\_{p,t}E_t\\ as well, and differentiating it requires
\\\partial^3 e_t/\partial u^3\\.

The derivative is only ever contracted. The full third derivative is an
\\m^3\\ array per observation and is never formed: the criterion needs
it only along the one direction in which the mode moves, so the method
propagates a matrix per observation, which has the same size as the
curvature and therefore the same \\O(nm^2)\\. A caller that needs
several hyperparameters calls this once per direction, which is cheaper
than one \\O(nm^3)\\ assembly whenever the hyperparameters are fewer
than the unknowns.

Each order of differentiation of the predictor involves one more order
of the family, because the recursion is driven by a score read at the
predictor that it produces: the curvature's `M` is built from third
derivatives, and this method needs a fourth. `blocks` therefore returns
two quantities beyond the curvature's, at an observation, from the
current Jacobian row \\D\\ and the direction \\v\\:

\$\$\texttt{dcurv} = \sum_r \ell\_{ppr}V_r, \qquad \texttt{N} =
\sum\_{r,r'}\Big(\sum\_{r''}\ell\_{prr'r''}(V\_{r''}\cdot v)\Big)
V_r^\top V\_{r'},\$\$

where `dcurv` serves for both the derivative of the curvature along
\\v\\, which is \\\texttt{dcurv}\cdot v\\, and the two terms that arise
from differentiating \\M\\'s own \\V_p\\.

The base method returns zeros, which is correct for a term whose
predictor is a block of columns times coefficients, since the second
derivative of such a predictor in the coefficients is zero. A term that
bends the predictor has to implement the method itself, because the base
method would report zeros for it. A structural term signals an error
instead: it has a second derivative, so a zero third derivative would
misstate it.

## See also

[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
for the second order this differentiates,
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
for the first,
[`statmodels7::reml()`](https://statmodels7.github.io/statmodels7/reference/reml.html)
for the criterion whose gradient requires it.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20))
term <- term_build(linpar(~t), dd)
# an additive term bends no predictor, so every order above the first is
# zero, as the base method returns
out <- term_third(term, rep(0, 20), dd$y,
                  score = function(e, i) dd$y[i] - e,
                  curvature = function(e, i) -1, psi = list(),
                  g = rep(1, 20), seed = matrix(0, 20, 2),
                  blocks = function(e, i, D, act) NULL,
                  direction = c(1, 0))
all(out$curvature == 0)
#> [1] TRUE
```
