# Where a Term's Own Coefficients Begin

Returns the starting coefficients of a built term, one per column of its
block. The base method returns zero everywhere, which is an admissible
start for a term whose block is a fixed design, since that block does
not depend on the coefficients.

## Usage

``` r
term_coef_start(term, target = NULL, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).
  An unbuilt one signals an error through
  [`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

- target:

  Optional numeric vector, one value per observation: the response on
  the scale of the predictor, `NULL` by default. It is supplied only
  where `params_interpretation` records that the response reads the
  parameter directly, so a term in a scale's equation receives `NULL`.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A numeric vector of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
in the block's column order.

## Why a term supplies a start of its own

A term that recomputes its block from its coefficients
([`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md))
needs a start of its own, because zero is degenerate for it. In
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
the break-point is read off two coefficients as \\-g_k/\delta_k\\, so a
vector of zeros puts every break-point at the same clamped position, and
with more than one break-point the columns of the block coincide. In
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
the Jacobian column is \\-\gamma_k \mathbb{1}(x \> \psi_k)\\ and
vanishes identically at \\\gamma_k = 0\\.

Without a `target`, the break-point terms return the start that
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
computed: unit changes, and the break-points at the interior quantiles
of the covariate or at the positions that `psi` names. With a `target`
of the right length and no developed coefficient, the slope and the
changes are instead the least-squares coefficients of the target on an
intercept and the term's columns at the starting positions.
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
returns the starting values of its own parameters carried through their
links.

A coefficient of zero has a meaning that depends on the block that the
term builds, so this start belongs to the term, as
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
does for a structural term.

## The role of `target`

`target` is the response carried onto the scale of the predictor to
which the term contributes. A term cannot compute it, because the term
holds its formula and its charts but not the distribution, the link or
the equation, which the fitting layer holds.
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) uses
it to estimate its own parameters from the data, over a deterministic
grid on each free parameter's chart, and the break-point terms use it
for the least-squares coefficients described above. Every other term
ignores it. The argument is optional; without it each term returns the
start that does not depend on the response.

## See also

[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
for a structural term's own parameters,
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for the terms for which this function exists,
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
for the grid rule behind a break-point start.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(50, 0, 10)))

# A fixed design starts at zero.
term_coef_start(term_build(linpar(~ x), dd))
#> [1] 0 0

# A break-point term does not: with several break-points, zero would
# place all of them at the same position.
jb <- term_build(jump(x, npsi = 2), dd)
setNames(term_coef_start(jb), term_coef_names(jb))
#> jump.delta1 jump.delta2     jump.g1     jump.g2 
#>    1.000000    1.000000   -3.847738   -6.908287 

# seg() starts the change at one and the break-point inside the data.
sb <- term_build(seg(x, npsi = 1), dd)
setNames(term_coef_start(sb), term_coef_names(sb))
#>   seg.beta seg.gamma1   seg.psi1 
#>   0.000000   1.000000   5.629448 
range(dd$x)
#> [1] 0.1339033 9.9190609
```
