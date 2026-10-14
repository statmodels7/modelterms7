# Coefficients a Term Draws for Itself

Replaces, in a draw of the term's block coefficients, the coefficients
whose scale the term determines. Every other coefficient is left as it
arrived, so a caller drawing a truth passes what it has drawn and
receives a vector that it can use.

## Usage

``` r
term_coef_draw(term, coef, sd = 1, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- coef:

  The coefficients drawn so far, of length
  [`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
  in the block's column order.

- sd:

  The width of the draw, `1` by default, multiplying the derived width
  above.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A list with `coef`, the coefficients with the term's own replaced;
`index`, the positions that were replaced, empty where none were; and
`scale`, one width per replaced position. A caller that reports the
truth of a prior over those coordinates reads `scale`, since the values
there are no longer a draw from that prior.

## Why a term draws at all

A coefficient of a design column is drawn from a normal of the width
that its caller chose, and for almost every term this is sufficient: a
slope is measured in the response's units against the covariate's, and
the term does not change that. A break-point is a position on the
covariate's own axis, so a draw made without reading the covariate has
no reason to fall inside the data. A draw from a standard normal falls
outside the range of a covariate that is measured on another scale, and
the confinement then moves the position to the edge of the interval,
where one of the two segments holds a twentieth of the rows and no break
is visible. The break-point terms therefore draw their own positions on
the covariate's axis.
[`term_draw()`](https://statmodels7.github.io/modelterms7/reference/term_draw.md)
applies the same division to the own parameters of a structural term:
the term that defines a quantity draws it.

## The width

The width follows from the spacing of the positions that the draw
spreads from. Writing \\w\\ for the width of the confinement interval
(from the 5th to the 95th percentile of the covariate) and \\K\\ for the
number of break-points, those positions are the \\j/(K+1)\\ quantiles of
the covariate. For a covariate spread evenly over the interval,
consecutive positions are \\w/(K+1)\\ apart. Requiring three standard
deviations to stay within half of that spacing, so that two neighbors
cross only if both travel further than a draw usually goes, gives

\$\$\sigma = \frac{sd \cdot w}{6(K+1)},\$\$

which at one break-point is \\w/12\\ and, for a covariate spread evenly
over the interval, keeps the draw inside \\\[0.25w, 0.75w\]\\ of the
interval at three standard deviations. For another distribution of the
covariate the spacing is a rule of thumb, and the width is still
computed from \\w\\.

## Where the coefficients go

A continuous or smoothed construction holds the break-point directly, so
the drawn position is written where it stands. A discontinuous one does
not hold it: it reads it off \\\psi_k = -g_k/\delta_k\\, so the
coefficient written is \\-\delta_k \psi_k\\ at the change of level that
the rest of the draw supplies. That is why the whole vector is passed in
and not only the term's own slice, and why the position is not the ratio
of two independently drawn normals, which has no mean and may fall
anywhere. For
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
the position read from the written coefficient is near the drawn one and
not equal to it.

Where the break-point carries a development the drawn positions are
projected onto that development's own design by least squares, so any
subformula is served and the spread the design induces is the one a
group actually gets.

## See also

[`term_draw()`](https://statmodels7.github.io/modelterms7/reference/term_draw.md)
for a structural term's own parameters,
[`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
for the point this spreads from,
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for the positions that it produces.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))

# For an ordinary term the draw is returned unchanged.
lb <- term_build(linpar(~ x), dd)
identical(term_coef_draw(lb, c(1, 2))$coef, c(1, 2))
#> [1] TRUE

# A break-point term places its own on the covariate's scale, where a
# draw from a standard normal would often fall outside the data.
sb <- term_build(seg(x, npsi = 1), dd)
d <- term_coef_draw(sb, stats::rnorm(term_npar(sb)))
setNames(d$coef, term_coef_names(sb))
#>    seg.beta  seg.gamma1    seg.psi1 
#> -0.62036668  0.04211587  5.16598496 
range(dd$x)
#> [1] 0.1307758 9.9268406
```
