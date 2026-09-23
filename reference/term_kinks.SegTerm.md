# Where a Break-Point Term's Objective Has a Kink

The positions, within the term's own coefficient vector, of the
coefficients that move a break-point sitting on an observation, for the
continuous construction only. There the truncated line \\(x-\psi)\_+\\
has a derivative in \\\psi\\ of \\-1\\ from one side and \\0\\ from the
other, so the objective is not differentiable in any coefficient that
moves \\\psi\\: the break-point's own when it carries no development,
and otherwise every coefficient of its development whose column is not
zero at the observations concerned. An observation counts as sitting on
the break-point within the band the block's own column reads its side
in.

A discontinuous construction answers `integer(0)`: its position is read
off a product of the unknowns and it is fitted by working fits rather
than by a line search on the objective. A smoothed term answers
`integer(0)`, having no kink, and so does a marginal one, which
integrates the break-point out.

## Arguments

- term:

  A built
  [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
  or
  [`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
  term.

- coef:

  The term's coefficients. `NULL` reads the stored ones.

- ...:

  Unused.

## Value

An integer vector of positions in the term's coefficients, possibly
empty.

## Examples

``` r
d <- data.frame(x = seq(0, 10, length.out = 41))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0)
b <- term_build(seg(x), d)
# the break-point placed on the observation at x = 6 is a kink in psi1
cf <- c(0.5, 2, 6)
term_kinks(b, cf)
#> [1] 3
# and between two observations it is not
term_kinks(b, c(0.5, 2, 6.1))
#> integer(0)
```
