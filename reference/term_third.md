# Third Derivatives of a Structural Term's Predictor, in One Direction

The second derivative of
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
differentiated once more along a single direction, and the derivative of
the term's Jacobian along that same direction. It is what the exact
gradient of a marginal criterion needs when a penalty covers the term's
own parameters.

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

  The weights the third derivative is contracted against, one per
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
direction the mode moves in. Where the predictor is \\X\beta\\ that
contraction is the family's third derivative against the design and
nothing else. Where a filter produces the predictor, \\K\\ carries
\\\sum_t w_t \ell\_{p,t}E_t\\ as well, and differentiating it asks for
\\\partial^3 e_t/\partial u^3\\.

**Only contracted.** The full third derivative is an \\m^3\\ array per
observation and is never formed: the criterion asks for it along the one
direction the mode moves in, so what is propagated is a matrix per
observation, the same size as the curvature and therefore the same
\\O(nm^2)\\. A caller wanting several hyperparameters calls this once
per direction, which is cheaper than one \\O(nm^3)\\ assembly whenever
the hyperparameters are fewer than the unknowns.

**What the model supplies.** Each order of differentiation of the
predictor pulls in one more order of the family, the recursion being
driven by a score read at the predictor it produces: the curvature's `M`
is built from third derivatives, and this needs a fourth. `blocks`
therefore returns two quantities beyond the curvature's, at an
observation with the current Jacobian row \\D\\ and its directional
derivative:

\$\$\texttt{dcurv} = \sum_r \ell\_{ppr}V_r, \qquad \texttt{N} =
\sum\_{r,r'}\Big(\sum\_{r''}\ell\_{prr'r''}(V\_{r''}\cdot v)\Big)
V_r^\top V\_{r'},\$\$

`dcurv` serving both the derivative of the curvature along \\v\\, which
is \\\texttt{dcurv}\cdot v\\, and the two terms differentiating \\M\\'s
own \\V_p\\.

The base method returns zeros. An additive term's second derivative is
already zero, so it is covered without writing anything, and a term
written later that does not implement this reports no third derivative
instead of a wrong one. A structural term refuses: it has a second
derivative, so zero there would be a false statement rather than a true
one.

## See also

[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
for the second order this differentiates,
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
for the first,
[`statmodels7::reml()`](https://statmodels7.github.io/statmodels7/reference/reml.html)
for the criterion whose gradient asks for it.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20))
term <- term_build(linpar(~t), dd)
# an additive term bends no predictor, so every order above the first is
# zero and the base method says so
out <- term_third(term, rep(0, 20), dd$y,
                  score = function(e, i) dd$y[i] - e,
                  curvature = function(e, i) -1, psi = list(),
                  g = rep(1, 20), seed = matrix(0, 20, 2),
                  blocks = function(e, i, D, act) NULL,
                  direction = c(1, 0))
all(out$curvature == 0)
#> [1] TRUE
```
