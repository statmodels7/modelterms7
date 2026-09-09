# A Break-Point Term's Own Coefficients, Drawn

The break-points drawn on the covariate's own axis and written back
through whichever route the construction reads them by, with every other
coefficient left as it arrived.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).
  An unbuilt one throws
  `"the term has not been built; call term_build(term, data) first."`.

- coef:

  The coefficients drawn so far, one per column of the block.

- sd:

  The width of the draw, `1` by default.

- ...:

  Unused.

## Value

A list with `coef`, `index` and `scale`, as
[`term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.md)
documents. `index` is every coordinate of every break-point and `scale`
the width each was drawn at, carrying the change of level as a factor
where the construction reads its position off one.

## Details

A break-point is a position on the covariate and nothing else in the
model is measured in those units, so it is the one coefficient here a
caller cannot draw. The positions are the covariate's own interior
quantiles plus a normal of the width
[`term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.md)
derives, confined to the interval between the 5th and the 95th
percentile.

Those quantiles and not the positions the term is carrying: they are
where
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
places a break-point nobody named, and they depend on the covariate
alone. A term handed over by a caller may sit somewhere else – a
simulation builds its specification against a placeholder response, and
a break-point chosen on a least-squares profile of noise lands against a
confinement limit as readily as anywhere – and a draw spreading from
there would inherit that.

A continuous or smoothed construction holds the position, so it is
written where it stands. A discontinuous one reads it off
\\-g_k/\delta_k\\, so what is written is \\-\delta_k \psi_k\\ at the
change of level the rest of the draw carries, which is
[`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md)'s
rule and is exact for a pure step, its read-off being memoryless. The
joint construction reads a quadratic that also carries the change of
slope and the increment left by the previous position, so there the
placement is the linear reading of it and the realized position is near
the drawn one rather than equal to it.

Where the break-point carries a development, the drawn positions are
projected onto that development's design, and the draw is made on its
coefficients scaled by the median row norm so that the spread arriving
at an observation is the derived width whatever the design looks like.

## See also

[`term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.md)
for the generic and the derivation of the width,
[`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md)
for placing the positions at named values.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(300, 0, 10)))
b <- term_build(seg(x, npsi = 1), d)

# Drawn on the covariate's scale, inside the confinement.
p <- replicate(50, {
  cf <- term_coef_draw(b, stats::rnorm(term_npar(b)))$coef
  seg_psi(term_refresh(b, cf))
})
range(p)
#> [1] 3.345020 6.668196
stats::quantile(d$x, c(0.05, 0.95), names = FALSE)
#> [1] 0.7697495 9.3469206
```
