# Build a Regime Term

Resolves the grouping and the ordering against the data and records them
in the blueprint. Nothing else is computed: the term has no design
block, and its chain and levels are already fixed by the constructor.

## Arguments

- term:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- data:

  A data frame carrying whatever `by` and `time` name.

- ...:

  Unused.

## Value

The term with `blueprint` filled.
[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
reads that property on this branch, so it is `TRUE` for the result.

## Details

`by` and `time` are evaluated in `data` with
[`baseenv()`](https://rdrr.io/r/base/environment.html) as the enclosure.
Each must give one value per row, and neither may be missing; both are
checked here, since this is the first point at which they meet data.

The blueprint holds `order`, the row indices of each group sorted by
time, and `n`, the observation count.
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
runs one forward recursion per group in that order, each from the
chain's stationary distribution, so a panel of independent series is
fitted by giving `by`.

Without `time` the rows are taken in the order they appear. That is a
real choice worth making deliberately: the recursion is about order, so
a data frame that is not already sorted gives a different model.

## See also

[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md),
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, id = rep(1:2, each = 20),
                 y = c(rnorm(20), rnorm(20, 3)))

# One group in row order.
b <- term_build(regime(2, time = t), dd)
names(b@blueprint)
#> [1] "order" "n"    
lengths(b@blueprint$order)
#>  1 
#> 40 

# Two independent series, each with its own recursion.
b2 <- term_build(regime(2, by = id, time = t), dd)
lengths(b2@blueprint$order)
#>  1  2 
#> 20 20 

# Both must give one value per row.
try(term_build(regime(2, time = c(1, 2)), dd))
#> Error : 'by' and 'time' must evaluate to one value per row.
```
