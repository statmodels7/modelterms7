# Reset a Break-Point Term's Scaling Schedule

The term with its scaling factors back at `c0`, the directions and the
step record cleared, and the break-points where they are. A restart
needs it: the schedule is a state of the iteration that only ever
tightens, so an iteration resumed from a converged fit inherits factors
at their floor and cannot travel. Measured, bootstrap restarts without
the reset returned the incumbent unchanged ten times out of ten.

## Usage

``` r
seg_reheat(term)
```

## Arguments

- term:

  A built break-point term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Value

The term, ready to iterate afresh from its current positions.

## See also

[`seg_step()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md),
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(100, 0, 10)))
b <- term_build(jump(x, psi = 4), dd)
b <- term_refresh(b, b@blueprint$coef)
identical(seg_reheat(b)@blueprint$cscale, b@blueprint$cscale * 0 + 0.05)
#> [1] TRUE
```
