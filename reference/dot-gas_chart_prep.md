# The Chart's Jets, Lifted to a Caller's Unknowns

Computes the level, the loadings, the autoregressive coefficients and
the starting level of a score-driven term, with their derivatives in a
caller's unknowns to the order that the number of directions sets: two
without a direction, the third contracted against one direction, and the
fourth contracted against two.

## Usage

``` r
.gas_chart_prep(zeta, p, q, links, zbase, act, dirs)
```

## Arguments

- zeta:

  The term's base parameters on the unconstrained scale.

- p, q:

  The score and autoregressive orders.

- links:

  The links, as
  [`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
  gives them.

- zbase:

  The positions of the base parameters among the unknowns.

- act:

  The unknowns over which the jets are formed.

- dirs:

  A list of zero, one or two directions over all the unknowns.

## Value

A list of the values and their lifted derivatives.

## Details

The chart's own derivatives come from
[`.gas_chart_derivs()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs.md)
and its siblings, in the term's base coordinates, and are lifted to the
unknowns by placing base coordinate \\j\\ at position `zbase[j]` of
`act`.
