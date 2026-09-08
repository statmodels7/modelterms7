# Drawing a Term's Own Parameters

A point drawn from the region a structural term's parameters plausibly
occupy: the values
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
gives, displaced by normal noise on the unconstrained scale
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

  The width of the draw, `1` by default, halved before it is applied for
  the reason measured above.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A named numeric vector on the unconstrained scale, of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
and named as
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
exactly as
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
returns.

## Details

Simulating from a model that carries state means choosing values for the
term's own parameters, and in quantity they cannot be chosen by hand: a
score-driven filter whose level is developed over a hundred groups
carries a hundred and one of them, and a caller made to name each one
names none and simulates a model with no dynamics and no heterogeneity
at all.

The rule is the start plus noise, and the start rather than zero because
the term already knows where the sensible region of its own charts is. A
score loading rides a log chart, so zero there is a loading of one and a
filter strong enough to destabilize its own recursion;
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
places it at 0.1 and the draw spreads from that. A constraint a chart
carries is then respected for free:
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
writes its levels as a first level and positive gaps, so any draw on the
unconstrained scale is ordered by construction, and a persistence on a
partial-autocorrelation chart is stationary at any coordinate whatever.

## The width

Half of `sd`, and the half is measured rather than chosen. Sixty draws
per family of a `gas(1, 1)` over sixty times, at widths 0.3, 0.5 and 1:
at a width of 1 two Poisson series of sixty are not finite and the level
of another spans 859, where at 0.5 none of a hundred and eighty fails
and the widest level spans 2.04. Read on the charts the same number is a
persistence between -0.27 and 0.88 at the fifth and ninety-fifth
percentiles, against -0.80 to 0.98 at a width of 1, which is the whole
stationary region and describes noise rather than a random effect.

## Coordinates a penalty covers

They are drawn here too, and a caller that can see the penalty is
expected to overwrite them. What a development's deviation is drawn from
is the prior its sub-term declares – a Gaussian random effect, a Laplace
one – and this method cannot see it: the prior's own scale is a
hyperparameter and choosing one is the caller's business.
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
says which coordinates those are, and penalties7's `penalty_draw()`
draws them.

## See also

[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
for the point it spreads from,
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the scale it is on,
[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md)
for the response drawn once these are chosen.

## Examples

``` r
set.seed(1)
g <- gas(p = 1, q = 1)

# The same names and the same scale as term_start(), drawn.
identical(names(term_draw(g)), term_params(g))
#> [1] TRUE

# Read on the parameters' own charts, a loading stays positive and a
# persistence stays stationary, whatever comes out.
lk <- term_links(g)
z <- term_draw(g)
c(alpha1 = linkfunctions7::linkinv(lk$alpha1, z[["alpha1"]]),
  pacf1 = linkfunctions7::linkinv(lk$pacf1, z[["pacf1"]]))
#>     alpha1      pacf1 
#>  0.1179103 -0.3886715 

# A wider draw is a wider model.
stats::sd(replicate(200, term_draw(g, sd = 2)[["omega"]]))
#> [1] 1.041263
```
