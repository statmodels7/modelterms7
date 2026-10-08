# The Parameters of a Marginal Break-Point Term

Numbered, one set per break-point, in a fixed order: the linear effect
`beta` where the kind carries one, then the prior's parameters for each
break-point, then the changes of slope `gamma1` ... and the changes of
level `delta1` ....

## Arguments

- term:

  A
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md),
  built or not.

- ...:

  Unused.

## Value

A character vector, of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

## Details

The prior's parameters are `psik.mean` and `psik.sd` under the default
Gaussian, the mean and the standard deviation of the positions of
break-point \\k\\ in the population. Where `random(distrib = )` named
another family the names are that family's own, its location fixed at
zero and `psi1.mean` carrying the position, so a Student t prior adds
`nu`.

Which of `gamma` and `delta` appear is the kind: `"seg"` has the changes
of slope, `"jump"` the changes of level, `"jseg"` both. Only `"seg"` and
`"jseg"` carry `beta`.

## See also

[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the chart each rides,
[`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md)
for what they mean.

## Examples

``` r
set.seed(1)
dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)

# A step term: the prior's location and scale, and the change of level.
term_params(term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd))
#> [1] "psi1.mean" "psi1.sd"   "delta1"   

# A continuous one adds the linear effect and the change of slope.
term_params(term_build(seg(x, psi ~ random(~ 1 | id), marginal = TRUE), dd))
#> [1] "beta"      "psi1.mean" "psi1.sd"   "gamma1"   
```
