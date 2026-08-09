# The Break-Points of a Segmented Term

The estimated positions of the break-points, one per break-point and per
level of `by`. For a continuous term these are coefficients; for a
discontinuous one they are read off two coefficients as \\-g/\kappa\\,
so this is the function that reports them either way.

## Usage

``` r
seg_psi(term, coef = NULL)
```

## Arguments

- term:

  A built
  [`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- coef:

  Optional coefficients; defaults to the ones the term was last
  refreshed at.

## Value

A numeric vector of break-point positions.

## See also

[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(100, 0, 10)))
seg_psi(term_build(seg(x, npsi = 2), dd))
#> [1] 3.823880 6.870228
```
