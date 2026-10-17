# The Break-Points of a Break-Point Term

Returns the estimated positions of the break-points. For a continuous or
a smoothed term these are coefficients; for a sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
they are read off two coefficients as \\-g_k/\delta_k\\, and for a sharp
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
off a quadratic in \\g_k\\, \\\delta_k\\, \\\gamma_k\\ and the previous
position, so this is the function that reports them under every
construction.

## Usage

``` r
seg_psi(term, coef = NULL)
```

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- coef:

  Optional coefficients; defaults to the coefficients at which the term
  was last refreshed.

## Value

A numeric vector with one entry per break-point, or, where a break-point
carries a development, a matrix with one row per observation and one
column per break-point.

## See also

[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(100, 0, 10)))
seg_psi(term_build(seg(x, npsi = 2), dd))
#> [1] 3.823880 6.870228
```
