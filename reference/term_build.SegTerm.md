# Build a Break-Point Term

Evaluates the covariate, resolves any development of the term's own
coefficients, places the break-points at their starting positions and
builds the working block there. It also settles the confinement limits,
the scaling schedule and, where the term is smoothed, the transition
width.

## Arguments

- term:

  A
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md),
  built or not.

- data:

  A data frame carrying the covariate and whatever the subformulas name.

- ...:

  Unused.

## Value

The term with `X`, `coef_names`, `blueprint` and, where a subformula
brought one, `penalty` filled.

## What is settled here

The covariate must vary. The break-points are confined to the interval
between the 5th and the 95th percentile of it: outside that the
indicator is constant, the truncated line and that constant are linearly
dependent, and the block is exactly singular. A run ending against the
limit has found no break-point, and
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
then reports the limit itself.

The starting positions are `psi` where it was given and evenly spaced
quantiles of the covariate otherwise. Zero is degenerate here, so
[`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
returns these rather than a vector of zeros.

Each subformula's right-hand side goes through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
and its terms are built, so their blueprints are recorded and reapplied
at prediction. What a development may carry depends on `kind`: the
continuous construction takes any of them, and the discontinuous ones
only a scalar change of level or a design that partitions the rows,
their read-off being a product of two unknowns.

## The smoothed form

With `smoothed` the step and the hinge are replaced by a
[`penalties7::abs_smoother()`](https://statmodels7.github.io/penalties7/reference/abs_smoother.html)'s
versions, and the transition width is resolved here from the covariate's
spacing, within groups where a development of the break-point supplies a
partition. The width is checked against the derived floor
\\\sqrt{\epsilon}D\\ and reported by
[`print()`](https://rdrr.io/r/base/print.html), a number that changes
what the term means having to be legible.

## See also

[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md);
[`term_refresh.SegTerm()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.SegTerm.md)
for the block as the break-points move;
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
for a better starting position than the quantiles.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(120, 0, 10)), id = factor(rep(1:4, each = 30)))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)

b <- term_build(seg(x, npsi = 1), d)
term_coef_names(b)
#> [1] "seg.beta"   "seg.gamma1" "seg.psi1"  
seg_psi(b)                     # the starting quantile
#> [1] 4.803127
setNames(term_coef_start(b), term_coef_names(b))
#>   seg.beta seg.gamma1   seg.psi1 
#>   0.000000   1.000000   4.803127 

# A developed break-point expands into its own design's coefficients.
bd <- term_build(seg(x, psi ~ id), d)
term_coef_names(bd)
#> [1] "seg.beta"             "seg.gamma1"           "seg.psi1.(Intercept)"
#> [4] "seg.psi1.id2"         "seg.psi1.id3"         "seg.psi1.id4"        
lapply(term_components(bd), function(z) z$index)
#> $beta
#> [1] 1
#> 
#> $gamma1
#> [1] 2
#> 
#> $psi1
#> [1] 3 4 5 6
#> 
```
