# Where a Nonlinear Term's Objective Has a Kink

The coefficients of the term that move a break-point sitting on an
observation, where a parameter is developed over a sharp
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md):
the answer of that sub-term's own
[`term_kinks()`](https://statmodels7.github.io/modelterms7/reference/term_kinks.md)
carried to the term's coefficients. Every other nonlinear term answers
`integer(0)`.

## Arguments

- term:

  A built
  [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).

- coef:

  The term's coefficients. `NULL` reads the stored ones.

- ...:

  Unused.

## Value

An integer vector of positions in the term's coefficients, possibly
empty.

## Examples

``` r
d <- data.frame(t = seq(0, 10, length.out = 41), x = rep(1:2, length.out = 41))
b <- term_build(nl(~ a * exp(-r * x), r ~ seg(t),
                   start = list(a = 2, r = 0.5)), d)
term_kinks(b)
#> [1] 5
```
