# Where a Term's Objective Has a Kink

The positions, among the term's own coefficients, of those in which the
term's contribution to the objective is not differentiable at the
coefficients given. A fitting layer reads them where a line search
rejects every step along a Newton or Gauss-Newton direction: holding
those coordinates at their values and solving for the rest is the step
at a kink, as holding a coordinate at a bound is the step at a boundary.

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

  Passed to methods. The shipped methods do not use it.

## Value

An integer vector of positions in the term's coefficients, possibly
empty.

## Details

A term whose contribution is a smooth function of its coefficients has
no kink, and that is the base method's answer. Two shipped terms answer
otherwise: the continuous break-point construction,
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
whose truncated line \\(x-\psi)\_+\\ is not differentiable in \\\psi\\
at an observation, and
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
which reports the kinks of a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
sub-term that develops one of its parameters. A smoothed
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and the discontinuous
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
answer `integer(0)`. The kinks of a penalty are not reported here;
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
reports them.

In a penalized fit of `seg(x, psi ~ random(~ 1 | id))`, the mode can put
one group's break-point exactly on an observation, where the one-sided
slopes differ. A scoring step on all the coefficients can then be
rejected at every length. Holding the coordinates that move that
break-point and solving for the rest gives a step that reaches the mode.

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
