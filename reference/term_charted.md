# Charted Coordinates of a Term

Names the free values of a term that are coordinates of a chart onto a
bounded set, so that a coordinate running to infinity reaches the
boundary of the quantity that it describes. A statistic read at such a
point requires a check that the free scale cannot provide: the
derivative on the free scale equals the derivative on the bounded scale
multiplied by the derivative of the chart, and the derivative of the
chart tends to zero at the edge, so the derivative on the free scale
vanishes there whatever the derivative on the bounded scale is.

## Usage

``` r
term_charted(term, ...)
```

## Arguments

- term:

  A model term. For
  [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) it
  must have been built.

- ...:

  Unused.

## Value

A character vector for a structural term, an integer vector for
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
named by parameter (`integer(0)` when the term has no charted
parameter), and `character(0)` otherwise.

## Details

For a structural term the result is a character vector of names from
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md):
every parameter whose link in
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
is not the identity, and for
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
also the additive log-ratios of the transition matrix, whose chart is
the one that
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html)
provides and which
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
therefore reports as the identity. For
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) it
is an integer vector of positions in the term's block, named by
parameter: the scalar parameters (those with no subformula) whose link
is not the identity. For every other term the result is an empty vector.

## See also

[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the charts,
[`statmodels7::statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.html)
for the check that reads this.

## Examples

``` r
term_charted(regime(2))
#> [1] "gap2"   "alr1.1" "alr2.1"
term_charted(gas(p = 1, q = 1))
#> [1] "kappa1" "pacf1" 

# A scalar parameter of nl() on a log link is charted; the position is
# its column in the term's block.
nb <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 1),
                    links = list(r = linkfunctions7::log_link())),
                 data.frame(x = seq(0.1, 2, length.out = 10)))
term_charted(nb)
#> r 
#> 2 
```
