# Where a Term's Own Coefficients Begin

The coefficients a built term asks to be started at, one per column of
its block. The base method returns zero everywhere. A term whose block
is a fixed design wants exactly that: the objective is convex in those
coordinates and the fit reaches the same optimum from anywhere.

## Usage

``` r
term_coef_start(term, target = NULL, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).
  An unbuilt one throws through
  [`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

- target:

  Optional numeric vector, one value per observation: the response on
  the scale of the predictor, `NULL` by default. Supplied only where
  `params_interpretation` says the response reads the parameter
  directly, so a term in a scale's equation is handed nothing.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A numeric vector of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
in the block's column order.

## Why any term needs a start of its own

A term that recomputes its block from its coefficients
([`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md))
is the case this exists for, because zero there is degenerate. In
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
the break-point is read off two coefficients as \\-g_k/\delta_k\\, so a
vector of zeros puts every break-point at the same clamped position and
makes the block singular. In
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
the Jacobian column is \\-\gamma_k \mathbb{1}(x \> \psi_k)\\ and
vanishes identically at \\\gamma_k = 0\\.

Those terms return the start
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
computed: unit changes, and the break-points at the interior quantiles
of the covariate or at the positions `psi` names.
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
returns the starting values of its own parameters carried through their
links.

Only the term knows what a coefficient of zero means for the block it
builds, which is the same reason
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
belongs to a structural term.

## What `target` is for

`target` is the response carried onto the scale of the predictor the
term contributes to. It is the one thing a term cannot work out for
itself: the term knows its formula and its charts, and the fitting layer
knows the distribution, the link and the equation.
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) uses
it to estimate its own parameters from the data, over a deterministic
grid on each free parameter's chart; every other term ignores it. It is
optional, so the default is what every term did before it existed.

## See also

[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
for a structural term's own parameters,
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for the terms this exists for,
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
for the grid rule behind a break-point start.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(50, 0, 10)))

# A fixed design starts at zero.
term_coef_start(term_build(linpar(~ x), dd))
#> [1] 0 0

# A break-point term does not: its block would be singular there.
jb <- term_build(jump(x), dd)
setNames(term_coef_start(jb), term_coef_names(jb))
#> jump.delta1     jump.g1 
#>    1.000000   -5.629448 

# seg() starts the change at one and the break-point inside the data.
sb <- term_build(seg(x, npsi = 1), dd)
setNames(term_coef_start(sb), term_coef_names(sb))
#>   seg.beta seg.gamma1   seg.psi1 
#>   0.000000   1.000000   5.629448 
range(dd$x)
#> [1] 0.1339033 9.9190609
```
