# The Coefficients a Term Holds

The positions, among the term's own coefficients, of those the term
holds at their values: the fitting layer keeps them out of the solve,
out of the information it inverts and out of the reading of the mode. A
term whose coefficients are all estimated by the derivative of the
objective holds none, which is the base method's answer.

## Usage

``` r
term_held(term, ...)
```

## Arguments

- term:

  A term.

- ...:

  Passed to methods.

## Value

An integer vector of positions in the term's coefficients, possibly
empty.

## Details

The one shipped term that answers otherwise is a
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term put in its held state by
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md).
There the break-point is fixed at a minimum of the profile objective,
which is constant between consecutive observations and so has no
derivative to estimate the position with, and the coefficient slot that
carried it in the working construction becomes a column of zeros.

## See also

[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md),
[`term_kinks()`](https://statmodels7.github.io/modelterms7/reference/term_kinks.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(100, 0, 10)))
d$y <- 1 + 2 * (d$x > 6) + rnorm(100, sd = 0.3)
b <- term_build(jump(x, psi = 6), d)
term_held(b)
#> integer(0)
term_held(seg_hold(b))
#> [1] 2
```
