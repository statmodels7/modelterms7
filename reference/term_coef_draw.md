# Coefficients a Term Draws for Itself

A simulation's draw of the term's block coefficients, with the ones
whose meaning only the term knows replaced. Every other coefficient is
left as it arrived, so a caller drawing a truth passes what it has drawn
and takes back a vector it can use.

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

  Passed to methods. No shipped method reads anything here.

## Value

A list with `coef`, the coefficients with the term's own replaced;
`index`, the positions it replaced, empty where it replaced none; and
`scale`, one width per replaced position. A caller reporting the truth
of a prior over those coordinates reads `scale`, the values there no
longer being that prior's draw.

## Why a term draws at all

A coefficient of a design column is drawn from a normal of the width its
caller chose, and for almost every term that is the whole story: a slope
is measured in the response's units against the covariate's, and the
term knows neither. A break-point is not such a quantity. It is a
position on the covariate's own axis, so a draw made without reading the
covariate lands outside the data most of the time, and the confinement
then pins it to the interval's edge, where one of the two segments holds
a twentieth of the rows and no break is visible.

Measured on a panel of fifty groups with the covariate uniform on \\(0,
1)\\: drawn from a standard normal the break-points run from \\-2.02\\
to \\1.65\\, and thirty-eight of the fifty end pinned against a
confinement limit with twelve strictly interior. This is the rule
[`term_draw()`](https://statmodels7.github.io/modelterms7/reference/term_draw.md)
states for a structural term's own parameters, read one level down:
whoever knows what a quantity means draws it.

## The width

Derived rather than chosen. Writing \\w\\ for the width of the
confinement interval and \\K\\ for the number of break-points, the
starting positions sit at the \\j/(K+1)\\ quantiles of it, so
consecutive ones are \\w/(K+1)\\ apart. Requiring three standard
deviations to stay within half of that spacing, so that two neighbours
cross only if both travel further than a draw usually goes, gives

\$\$\sigma = \frac{sd \cdot w}{6(K+1)},\$\$

which at one break-point is \\w/12\\ and keeps the draw inside
\\\[0.25w, 0.75w\]\\ of the interval at three standard deviations.

## Where the coefficients go

A continuous or smoothed construction holds the break-point directly, so
the drawn position is written where it stands. A discontinuous one does
not hold it at all: it reads it off \\\psi_k = -g_k/\delta_k\\, so the
coefficient written is \\-\delta_k \psi_k\\ at the change of level the
rest of the draw supplies. That is why the whole vector is passed in and
not only the term's own slice, and why a ratio of two independently
drawn normals – which has no mean and puts the position anywhere – is
not what comes out.

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
for the positions it produces.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))

# An ordinary term has no opinion and hands the draw back.
lb <- term_build(linpar(~ x), dd)
identical(term_coef_draw(lb, c(1, 2))$coef, c(1, 2))
#> [1] TRUE

# A break-point term places its own on the covariate's scale, where a
# standard normal would have put it outside the data altogether.
sb <- term_build(seg(x, npsi = 1), dd)
d <- term_coef_draw(sb, stats::rnorm(term_npar(sb)))
setNames(d$coef, term_coef_names(sb))
#>    seg.beta  seg.gamma1    seg.psi1 
#> -0.62036668  0.04211587  5.16598496 
range(dd$x)
#> [1] 0.1307758 9.9268406
```
