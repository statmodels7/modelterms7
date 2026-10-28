# The Constraints That a Break-Point on an Observation Places

Returns, for each break-point of a sharp
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
that sits on an observation at the given coefficients, the row of the
linear constraint that keeps it there. The row is the unit vector of its
position where the position is not developed, and the development's row
\\z_j\\ at the position's coefficients where it is, so that the kink
lies on the hyperplane \\z_j^\top\kappa = x_j\\.

## Arguments

- term:

  A built
  [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
  term.

- coef:

  The term's coefficients. `NULL` reads the stored ones.

- ...:

  Unused.

## Value

A numeric matrix with one column per coefficient of the term and one row
per distinct constraint. If no break-point sits on an observation, the
matrix has no rows. A smoothed
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
has no rows.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(60, 0, 10), w = runif(60))
b <- term_build(seg(x, psi ~ w), d)
ix <- b@build_info$index$psi1
cf <- b@build_info$coef
# the break-point of the first observation is placed on that observation
cf[ix[1L]] <- d$x[1L] - cf[ix[2L]] * d$w[1L]
term_kink_constraints(b, cf)
#>      [,1] [,2] [,3]      [,4]
#> [1,]    0    0    1 0.9128759
```
