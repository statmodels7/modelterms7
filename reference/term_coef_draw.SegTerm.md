# Drawing the Coefficients of a Break-Point Term

Draws the break-points on the covariate's own axis and writes them back
through the route by which the construction reads them, whichever it is,
and every other coefficient is left as it arrived. A term held by
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md)
draws nothing: the coefficients are returned unchanged, with an empty
`index` and an empty `scale`, because the positions are held.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).
  An unbuilt one signals the error
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
the width at which each was drawn, carrying the change of level as a
factor where the construction reads its position off one. For a held
term `coef` is the argument unchanged, and `index` and `scale` are
empty.

## Details

A break-point is a position on the covariate and is the only coefficient
of the term measured in those units, so a draw from a standard normal
does not apply to it. The positions are the covariate's own interior
quantiles plus a normal of the width that
[`term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.md)
derives, confined to the interval between the 5th and the 95th
percentile.

The draw spreads from these quantiles and not from the positions that
the term currently holds.
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
places a break-point that has not been named at these quantiles, and
they depend only on the covariate. A term handed over by a caller may
hold positions chosen elsewhere: a simulation builds its specification
against a placeholder response, and a break-point chosen on a
least-squares profile of noise lands against a confinement limit as
readily as anywhere, so a draw spreading from there would inherit that
position.

A continuous or smoothed construction holds the position, so it is
written where it stands. A sharp discontinuous one reads the position
off \\-g_k/\delta_k\\, so the method writes \\-\delta_k \psi_k\\, taking
\\\delta_k\\ from the rest of the draw. This is the rule of
[`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md),
and it is exact for a pure step because the read-off of a pure step is
memoryless. The joint construction reads a quadratic that also carries
the change of slope and the increment left by the previous position, so
there the placement is the linear reading of it and the realized
position is near the drawn one and not equal to it.

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

# The positions are drawn on the covariate's scale, inside the confinement.
p <- replicate(50, {
  cf <- term_coef_draw(b, stats::rnorm(term_npar(b)))$coef
  seg_psi(term_refresh(b, cf))
})
range(p)
#> [1] 3.345020 6.668196
stats::quantile(d$x, c(0.05, 0.95), names = FALSE)
#> [1] 0.7697495 9.3469206
```
