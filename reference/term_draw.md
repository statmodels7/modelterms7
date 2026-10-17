# Drawing a Term's Own Parameters

Draws a point from the region that a structural term's parameters
plausibly occupy. The point consists of the values that
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
gives, displaced by normal noise on the unconstrained scale that
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
defines.

## Usage

``` r
term_draw(term, sd = 1, ...)
```

## Arguments

- term:

  A built structural term.

- sd:

  The width of the draw, `1` by default. The standard deviation applied
  is `sd / 2`, for the reason given in the details.

- ...:

  Passed to methods. The shipped methods do not use these arguments.

## Value

A named numeric vector on the unconstrained scale, of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
and named as
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
exactly as
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
returns.

## Details

Simulating from a model that carries state requires values for the
term's own parameters. A score-driven filter whose level is developed
over a hundred groups has a hundred level coordinates in addition to its
loading and persistence coordinates, and choosing each of them by hand
is impractical.

The rule is the start plus noise. The start is used instead of zero
because it lies in the region of the term's own charts where the term is
well behaved. A score loading rides a log chart, so zero there is a
loading of one and a filter strong enough to destabilize its own
recursion;
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
places it at 0.1 and the draw spreads from that. A constraint that a
chart carries is respected by every draw:
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
writes its levels as a first level and positive gaps, so any draw on the
unconstrained scale is ordered by construction, and a persistence on a
partial-autocorrelation chart is stationary at any coordinate.

## The width

The standard deviation of the draw is `sd / 2`. At the full width of
`sd` a score-driven filter with a Poisson response can produce
non-finite series, and a persistence on the partial-autocorrelation
chart covers nearly the whole stationary region, which describes noise
and not a random effect.

## Coordinates covered by a penalty

These coordinates are drawn as well, and a caller that has access to the
penalty is expected to overwrite them. The deviation of a development
follows the prior that its sub-term declares (a Gaussian random effect
or a Laplace one). This method does not receive that prior, because its
scale is a hyperparameter that the caller chooses.
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
identifies those coordinates, and penalties7's `penalty_draw()` draws
them.

## See also

[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
for the point from which the draw spreads,
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the scale on which the draw is expressed,
[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md)
for the response drawn once these are chosen.

## Examples

``` r
set.seed(1)
g <- gas(p = 1, q = 1)

# The draw has the same names and the same scale as term_start().
identical(names(term_draw(g)), term_params(g))
#> [1] TRUE

# Read on the parameters' own charts, a loading stays positive and a
# persistence stays stationary, whatever comes out.
lk <- term_links(g)
z <- term_draw(g)
c(kappa1 = linkfunctions7::linkinv(lk$kappa1, z[["kappa1"]]),
  pacf1 = linkfunctions7::linkinv(lk$pacf1, z[["pacf1"]]))
#>     kappa1      pacf1 
#>  0.1179103 -0.3886715 

# A larger sd widens the draws, and with them the range of models.
stats::sd(replicate(200, term_draw(g, sd = 2)[["omega"]]))
#> [1] 1.041263
```
