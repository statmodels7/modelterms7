# Confinement Intervals of the Break-Points of a Term

Returns the interval in which each break-point of a built
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term can lie: the interval between the 5th and the 95th percentile of
the covariate for every break-point, or the narrower interval that
[`seg_pin()`](https://statmodels7.github.io/modelterms7/reference/seg_pin.md)
gave to a break-point that is free beside a held one.

## Usage

``` r
seg_limits(term)
```

## Arguments

- term:

  A built break-point term.

## Value

A numeric matrix with one row for each break-point and the columns
`lower` and `upper`.

## See also

[`seg_pin()`](https://statmodels7.github.io/modelterms7/reference/seg_pin.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(100, 0, 10)))
seg_limits(term_build(seg(x, npsi = 2), d))
#>          lower    upper
#> psi1 0.8356852 9.084412
#> psi2 0.8356852 9.084412
```
