# The Constraints a Break-Point on an Observation Places

For each break-point of a sharp
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
that sits on an observation at the coefficients given, the row of the
linear constraint that keeps it there: the unit vector of its position
where the position is not developed, and the development's row \\z_j\\
at the position's coefficients where it is, the kink lying on the
hyperplane \\z_j^\top\kappa = x_j\\.

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
per distinct constraint, with no rows where no break-point sits on an
observation. A smoothed
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
has no rows.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(60, 0, 10), w = runif(60))
b <- term_build(seg(x, psi ~ w), d)
ix <- b@blueprint$index$psi1
cf <- b@blueprint$coef
# the break-point of the first observation placed on it
cf[ix[1L]] <- d$x[1L] - cf[ix[2L]] * d$w[1L]
term_kink_constraints(b, cf)
#>      [,1] [,2] [,3]      [,4]
#> [1,]    0    0    1 0.9128759
```
