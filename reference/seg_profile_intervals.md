# The Profile of One Break-Point Over Every Interval

Returns the least-squares profile of a sharp
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term in one of its break-points, the others held at their current
positions. The profile is the residual sum of squares in every interval
between consecutive distinct values of the covariate inside the
confinement limits, evaluated at the midpoint of the interval for a
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
and at the minimizing position inside the interval for a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md).
This is the quantity that
[`seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.md)
minimizes. It is returned whole so that a caller whose objective the
profile only approximates can take several candidates from it.

## Usage

``` r
seg_profile_intervals(term, y, k = 1L, weights = NULL)
```

## Arguments

- term:

  A built break-point term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- y:

  A numeric vector of finite values, one per observation of the build
  data: the response, net of whatever the caller wants held.

- k:

  Which break-point, an integer.

- weights:

  Optional non-negative weights, one per observation; the profile is
  then weighted least squares.

## Value

A data frame with columns `psi`, the position in each interval, and
`rss`, the profile there (`Inf` where the design is singular). The data
frame has no rows when no interval lies inside the confinement limits.

## See also

[`seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))
dd$y <- 1 + 1.5 * (dd$x > 6) + rnorm(200, sd = 0.4)
b <- term_build(jump(x, psi = 4), dd)
pr <- seg_profile_intervals(b, dd$y)
pr$psi[which.min(pr$rss)]
#> [1] 6.005535
```
