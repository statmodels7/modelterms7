# Polish a Break-Point Term's Positions on the Exact Profile

Coordinate descent over the break-point positions on the exact profile:
one position at a time is swept over a grid with the others held, the
least-squares fit at fixed positions being an ordinary linear model, and
the sweeps repeat until no position moves. The term comes back relocated
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

  A numeric vector, one value per observation of the build data: the
  response, net of whatever the caller wants held.

- k:

  How many grid points per sweep. Defaults to 50.

- sweeps:

  At most how many passes over the positions. Defaults to 10; the
  descent usually stops moving in two or three.

- weights:

  Optional non-negative weights, one per observation; the profile is
  then weighted least squares.

## Value

The term at the polished positions.

## Details

This is
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)'s
device made conditional, and it exists for the failure a start cannot
fix: with several break-points the iteration can capture two of them on
one feature and press the third against its confinement limit, and
neither a bootstrap excursion – whose perturbation is of order
\\1/\sqrt{n}\\, nor a random relocation escapes it, the basin of the
right configuration being narrow. A grid sweep of one position with the
others held walks straight to the missing feature, because the profile
at fixed positions is exact and unimodal around it.

The profile is least squares of `y` on the term's own columns at fixed
positions, so it is exact for a gaussian response and a starting rule
for any other, the argument
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
already makes. A term whose per-break-point coefficients carry a
development is rejected, its positions being one per observation.

With `weights`, the profile is weighted least squares. A restarting loop
sweeps the profile of a bootstrap resample that way, the multinomial
counts of sampling the rows with replacement as weights – which moves
the profile's optima the way refitting the resample would, at the cost
of grid-many linear fits instead of a whole model fit: the non-convexity
of these models lives entirely in the positions, so exploring the
positions is the whole of what an excursion is for.

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
