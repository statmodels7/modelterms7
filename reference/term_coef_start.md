# Where a Term's Own Coefficients Begin

The coefficients a built term asks to be started at, one per column of
its block. The base method returns zero everywhere, which is what a term
whose block is a fixed design wants: the fit reaches the same optimum
from anywhere, the objective being convex in those coordinates.

## Usage

``` r
term_coef_start(term, target = NULL, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- target:

  Optional numeric vector, one value per observation: the response on
  the scale of the predictor.
  [`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md) uses
  it to estimate its own parameters; every other term ignores it.

- ...:

  Passed to methods.

## Value

A numeric vector of length
[`term_npar`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

## Details

A term that recomputes its block from its coefficients
([`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md))
is the case this exists for, because zero is not a neutral point there
but a degenerate one. In
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)
the break-point is read off two coefficients as \\-g_k/\delta_k\\, so a
vector of zeros puts every break-point at the same clamped position and
makes the block singular; in
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md) the
Jacobian column is \\-\gamma_k\\\mathbb{1}(x \> \psi_k)\\ and vanishes
identically. Those terms return the start
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
computed – unit changes and the break-points at the interior quantiles
of the covariate, or the positions `psi` names – and
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md)
returns the starting values of its own parameters carried through their
links.

The value belongs to the term for the reason
[`term_start`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
records for a structural one: only the term knows what a coefficient of
zero means for the block it builds.

`target` is the response carried onto the scale of the predictor the
term contributes to, which is what a term needs to estimate parameters
of its own from the data and is the one thing it cannot work out for
itself: the term knows its formula and its charts, the fitting layer
knows the distribution, the link and the equation. It is optional, and a
term that has no use for it ignores it, so the default is the behaviour
every term had before it existed.

## See also

[`term_start`](https://statmodels7.github.io/modelterms7/reference/term_start.md),
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)

## Examples

``` r
dd <- data.frame(x = sort(runif(50, 0, 10)))
term_coef_start(term_build(linpar(~x), dd))
#> [1] 0 0
term_coef_start(term_build(jump(x), dd))
#> [1]  1.000000 -4.774248
```
