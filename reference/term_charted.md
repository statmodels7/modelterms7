# The Coordinates of a Term That Ride a Chart With an Edge

Names the free values of a term that are coordinates of a chart mapping
onto a bounded set, where a coordinate running to infinity reaches the
boundary of the quantity it describes. A statistic read at such a point
needs a check the free scale cannot give: the derivative on the free
scale is the derivative on the bounded scale times the chart's own
derivative, and the second vanishes at the edge whatever the first says.

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

A character vector for a structural term, a named integer vector for
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md), and
`character(0)` otherwise.

## Details

For a structural term the answer is a character vector of names from
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md):
every parameter whose link in
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
is not the identity, and for
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
also the additive log-ratios of the transition matrix, whose chart is
the one
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html)
provides and which
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
therefore reports as the identity. For
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) it
is an integer vector of positions in the term's block, named by
parameter: the scalar parameters (those with no subformula) whose link
is not the identity. Every other term answers with an empty vector.

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
```
