# Build a Break-Point Term

Evaluates the covariate, resolves any development of the term's own
coefficients, places the break-points at their starting positions and
builds the working block there. It also sets the confinement limits, the
scaling schedule and, where the term is smoothed, the transition width.

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

The term with `X`, `coef_names` and `blueprint` filled. The slot
`penalty` is `NULL`; the penalties of the developments are recorded in
the blueprint and listed by
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md).

## The quantities set at build

The covariate must vary. The break-points are confined to the interval
between the 5th and the 95th percentile of it, so that both sides of a
break-point carry observations (a break-point beyond the range of the
data would make the indicator constant, the truncated line and that
constant linearly dependent, and the block exactly singular). A run
ending against the limit has found no break-point, and
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
then reports the limit itself.

The starting positions are `psi` where it was given and evenly spaced
quantiles of the covariate otherwise. A start of zero is degenerate for
this term, so
[`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
returns these positions and not a vector of zeros.

Each subformula's right-hand side goes through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
and its terms are built, so their blueprints are recorded and reapplied
at prediction. The developments that a term may carry depend on `kind`.
The continuous construction takes any of them. A sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
takes a development of the break-point with a scalar change of level on
any design, or one development shared by the change of level and the
break-point on a design that partitions the rows; a sharp
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
takes only the shared partition development, because the read-off is a
product of two unknowns. A smoothed term takes any development. A design
that carries the intercept twice, the subformula's own and a constant
column of a sub-term, as in `gamma ~ lasso(~1)`, is rejected, and
`gamma ~ 0 + lasso(~1)` is the spelling of an intercept-only penalized
development.

## The smoothed form

With `smoothed` the step and the hinge are replaced by the versions of a
[`numericals7::abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.html),
and the transition width is resolved here from the covariate's spacing,
within groups where a development of the break-point supplies a
partition. The width is checked against the floor \\\sqrt{\epsilon}D\\,
with \\D\\ the range of the covariate, carried onto the smoother's own
width scale by
[`numericals7::smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.html),
and a width below it is rejected. The width is reported by
[`print()`](https://rdrr.io/r/base/print.html).

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
