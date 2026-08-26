# The Exact Profile of a Break-Point Term at Its Positions

The residual sum of squares of `y` on the term's own columns at its
current positions, the number
[`seg_polish()`](https://statmodels7.github.io/modelterms7/reference/seg_polish.md)
descends on, read at one point. A restarting loop screens proposals with
it: two configurations of positions are compared by their profiles at
the cost of two linear fits, where refitting the model to compare them
costs a whole fit each.

## Usage

``` r
seg_profile_rss(term, y, weights = NULL)
```

## Arguments

- term:

  A built break-point term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- y:

  A numeric vector, one value per observation of the build data: the
  response, net of whatever the caller wants held.

- weights:

  Optional non-negative weights, one per observation; the profile is
  then weighted least squares.

## Value

A single number; `Inf` where the fit fails.

## See also

[`seg_polish()`](https://statmodels7.github.io/modelterms7/reference/seg_polish.md),
[`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(100, 0, 10)))
dd$y <- 2 * (dd$x > 6) + rnorm(100, sd = 0.3)
b <- term_build(jump(x, psi = 4), dd)
seg_profile_rss(b, dd$y) > seg_profile_rss(seg_relocate(b, 6), dd$y)
#> [1] TRUE
```
