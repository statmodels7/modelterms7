# Polish a Break-Point Term's Positions on the Exact Profile

Performs coordinate descent over the break-point positions on the exact
profile: one position at a time is swept over a grid while the others
are held (the least-squares fit at fixed positions is an ordinary linear
model), and the sweeps repeat until a full sweep leaves every position
unchanged. The grid has `k` quantiles of the covariate, from the 2nd to
the 98th percentile and confined to the term's limits. The term comes
back relocated
([`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md))
at the best positions found.

## Usage

``` r
seg_polish(term, y, k = 50, sweeps = 10, weights = NULL)
```

## Arguments

- term:

  A built break-point term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- y:

  A numeric vector of finite values, one per observation of the build
  data: the response, net of whatever the caller wants to hold fixed.

- k:

  How many grid points per sweep. Defaults to 50.

- sweeps:

  At most how many passes over the positions. Defaults to 10. The
  descent stops once a full pass moves no position.

- weights:

  Optional non-negative weights, one per observation; the profile is
  then weighted least squares.

## Value

The term at the polished positions.

## Details

This is the device of
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
made conditional on the current positions, and it addresses the failure
that a start cannot fix: with several break-points the iteration can
capture two of them on one feature and press the third against its
confinement limit, and a bootstrap excursion, whose perturbation is of
order \\1/\sqrt{n}\\, and a random relocation both fail to escape it,
because the basin of the right configuration is narrow. A grid sweep of
one position, while the others are held, reaches the missing feature,
since the profile at fixed positions is exact.

The profile is least squares of `y` on an intercept, the linear column
where the term has one, and the term's own columns at fixed positions,
so it is exact for a gaussian response and a starting rule for any
other, the argument that
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
already makes. A term whose per-break-point coefficients carry a
development is rejected, because its positions are one per observation.

With `weights`, the profile is weighted least squares. A restarting loop
can sweep the profile of a bootstrap resample in this way, using the
multinomial counts of sampling the rows with replacement as weights.
This moves the optima of the profile as refitting the resample would, at
the cost of one linear fit per grid point instead of a whole model fit.

## See also

[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md),
[`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(300, 0, 10)))
dd$y <- 2 * (dd$x > 3) - 1.5 * (dd$x > 7) + rnorm(300, sd = 0.3)
b <- term_build(jump(x, npsi = 2, psi = c(1, 2)), dd)
seg_psi(seg_polish(b, dd$y))
#> [1] 2.937291 7.011637
```
