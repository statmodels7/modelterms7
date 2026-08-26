# Build a Marginal Break-Point Term

Evaluates the covariate and the grouping variable, works out the
interval structure the marginal likelihood is summed or integrated over,
and records both in the blueprint together with a data-based starting
point for the term's parameters.

## Arguments

- term:

  A
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- data:

  A data frame carrying the covariate and the grouping variable.

- ...:

  Unused.

## Value

The term with `blueprint` filled.
[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
reads that property on this branch, so it is `TRUE` for the result.

## Details

The covariate must vary; a constant one has no break-point to place and
throws. The grouping variable must give one value per row.

What the blueprint records is the group of each observation, the
covariate sorted within each group, and the labels of the groups. The
intervals the likelihood decomposes over are the gaps between a group's
ordered observations, so they follow from that ordering.

The start is computed here because zero is degenerate for this term:
with `delta = 0` the intervals are indistinguishable, every mass
derivative sums to the derivative of a constant, and the surface is
exactly flat in the prior's location and scale.
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
returns what this computed.

## See also

[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md),
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)

b <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
names(b@blueprint)
#> [1] "n"      "x"      "groups" "labels"
b@blueprint$labels
#> [1] "1" "2" "3"

# The start is data-based, zero being degenerate here.
round(term_start(b), 4)
#>     m1   tau1 delta1 
#> 4.5000 0.1542 0.0000 

# A constant covariate has no break-point to place.
flat <- transform(dd, x = 1)
try(term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), flat))
#> Error : the covariate of a break-point term must vary.
```
