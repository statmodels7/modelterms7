# Build a Marginal Break-Point Term

Evaluates the covariate and the grouping variable, orders the rows of
each group by the covariate, and records the result in the blueprint.
The interval structure over which the marginal likelihood is summed or
integrated follows from that ordering.

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

The covariate must vary: a constant covariate has no break-point to
place and is rejected. The grouping variable must give one value per row
and must not contain missing values.

The blueprint records the number of rows (`n`), the covariate in the
order of the data (`x`), for each group the row indices ordered by the
covariate (`groups`), and the labels of the groups (`labels`). The
intervals over which the likelihood decomposes are the gaps between a
group's ordered observations, so they follow from that ordering.

No starting values are computed here:
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
computes them when it is called. A start with every change at zero is
degenerate for this term. With `delta = 0` the intervals are
indistinguishable, every mass derivative sums to the derivative of a
constant, and the surface is exactly flat in the location and scale of
the prior. A fitting layer therefore supplies a target to
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md),
from which the changes are estimated.

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

# Without a target the start uses the covariate's quantiles and spread, and
# the change of level starts at zero.
round(term_start(b), 4)
#> psi1.mean   psi1.sd    delta1 
#>    4.5000    0.1542    0.0000 

# A constant covariate has no break-point to place.
flat <- transform(dd, x = 1)
try(term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), flat))
#> Error : the covariate of a break-point term must vary.
```
