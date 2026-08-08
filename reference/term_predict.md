# Design Block on New Data

Applies a built term's mapping to new data, reproducing the block the
term would have produced had the new rows been part of the original
data: factor levels, contrasts and any constants recorded in the
blueprint at build time are reused, never recomputed. New data carrying
a factor level unknown to the blueprint is rejected.

## Usage

``` r
term_predict(term, newdata, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- newdata:

  A data frame.

- ...:

  Passed to methods.

## Value

A numeric matrix with `nrow(newdata)` rows and one column per
coefficient.

## Examples

``` r
built <- term_build(linpar(~x), data.frame(x = 1:4))
term_predict(built, data.frame(x = c(0.5, 2.5)))
#>   (Intercept)   x
#> 1           1 0.5
#> 2           1 2.5
```
