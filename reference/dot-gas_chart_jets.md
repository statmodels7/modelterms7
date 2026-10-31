# The Chart's Jets in the Caller's Unknowns

Returns the jets of the level, of each score loading, of each
autoregressive coefficient and of the starting level, over the given
directions, in the caller's unknowns.

## Usage

``` r
.gas_chart_jets(zeta, p, q, links, zbase, m, dirs)
```

## Arguments

- zeta:

  The base parameters on the unconstrained scale.

- p, q:

  The orders.

- links:

  The links.

- zbase:

  The columns of the unknowns that hold the base coordinates.

- m:

  The number of unknowns.

- dirs:

  The directions, in the unknowns.

## Value

A list with `omega`, `a` (a list), `b` (a list) and `f0`, jets.

## Details

The chart's derivatives are those of
[`.gas_chart_contract()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_contract.md)
in the base coordinates, lifted to the unknowns through the columns
`zbase` that hold the base coordinates. The starting level \\f_0 =
\omega/(1 - S)\\, with \\S\\ the sum of the autoregressive coefficients,
is read off the fixed point \\f_0 = \omega + S f_0\\: along each subset
of the directions, taken by increasing size, its matrix is \\(1 -
S)^{-1}\\ times the level's plus the product rule of \\S f_0\\ without
the term that carries every derivative on \\f_0\\.
