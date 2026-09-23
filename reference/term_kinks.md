# Where a Term's Objective Has a Kink

The positions, among the term's own coefficients, of those in which the
model's objective is not differentiable at the coefficients given. A
fitting layer reads it where a line search rejects every step along a
Newton or Gauss-Newton direction: holding those coordinates at their
values and solving for the rest is the step at a kink, the same thing
holding a coordinate at a bound is at a boundary.

## Usage

``` r
term_kinks(term, coef = NULL, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The term's coefficients, `NULL` for the stored ones.

- ...:

  Passed to methods.

## Value

An integer vector of positions in the term's coefficients, possibly
empty.

## Details

A term whose contribution is a smooth function of its coefficients has
no kink, and that is the base method's answer. The one shipped term that
answers otherwise is the continuous break-point construction,
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
whose truncated line \\(x-\psi)\_+\\ is not differentiable in \\\psi\\
at an observation. Measured on `seg(x, psi ~ random(~1 | id))`, the
penalized mode puts one group's break-point exactly on an observation,
the one-sided slopes there being \\+6.0\\ and \\-1.6\\, and a scoring
step on all the coefficients is rejected at every length: the other
groups' deviations then stay where they were, with scores up to \\1.6\\,
and the objective 0.010 above the mode. Holding the two coordinates that
move that break-point reaches the mode in four iterations.

## See also

[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md),
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md).

## Examples

``` r
d <- data.frame(x = seq(0, 10, length.out = 41))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0)
term_kinks(term_build(linpar(~ x), d))
#> integer(0)
term_kinks(term_build(seg(x), d), c(0.5, 2, 6))
#> [1] 3
```
