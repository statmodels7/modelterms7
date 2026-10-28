# Build a Regime Term

Resolves the grouping and the ordering against the data and records them
in the `build_info`. The method does not compute anything else: the term
has no design block, its chain structure is fixed by the constructor,
and its levels are parameters supplied at each evaluation.

## Arguments

- term:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- data:

  A data frame carrying whatever `by` and `time` name.

- ...:

  Unused.

## Value

The term, in which `build_info` is filled.
[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
reads that property on this branch, so it is `TRUE` for the result.

## Details

`by` and `time` are evaluated in `data` with
[`baseenv()`](https://rdrr.io/r/base/environment.html) as the enclosure.
Each must give one non-missing value per row; both are checked here,
since this is the first point at which they meet data.

The `build_info` holds `order`, the row indices of each group sorted by
time, and `n`, the observation count.
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
runs one forward recursion per group in that order, each from the
chain's stationary distribution, so a panel of independent series is
fitted by giving `by`.

Without `time` the rows are taken in the order in which they appear. The
recursion depends on that order, so a data frame that is not already
sorted defines a different model. With `time`, the rows of a group are
sorted by it, and `time` only orders them: a gap between two times
counts as one step of the chain.

## See also

[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md),
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, id = rep(1:2, each = 20),
                 y = c(rnorm(20), rnorm(20, 3)))

# The rows form one group in row order.
b <- term_build(regime(2, time = t), dd)
names(b@build_info)
#> [1] "order" "n"    
lengths(b@build_info$order)
#>  1 
#> 40 

# The rows form two independent series, and each has its own recursion.
b2 <- term_build(regime(2, by = id, time = t), dd)
lengths(b2@build_info$order)
#>  1  2 
#> 20 20 

# Both must give one value per row.
try(term_build(regime(2, time = c(1, 2)), dd))
#> Error : 'by' and 'time' must evaluate to one value per row.
```
