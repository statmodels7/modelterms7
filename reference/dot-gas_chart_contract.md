# The Chart's Derivatives of Any Order up to Five, in a Set of Directions

Returns the chart's derivative of order `length(dirs) + 2` in the term's
base coordinates, contracted against the directions in `dirs` and
leaving two slots free: one matrix for the level, one per score loading
and one per autoregressive coefficient. With zero, one and two
directions it equals the second, third and fourth orders of
[`.gas_chart_derivs()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs.md),
[`.gas_chart_derivs3()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs3.md)
and
[`.gas_chart_derivs4()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs4.md).

## Usage

``` r
.gas_chart_contract(zeta, p, q, links, dirs)
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

- dirs:

  A list of zero to three directions in the coordinates of `zeta`.

## Value

A list with `omega`, `a` and `b`: a matrix, a list of `max(p, 1)`
matrices and a list of `q` matrices over the base coordinates.

## Details

The level and the loadings ride scalar links, so their derivative is a
single diagonal entry, the link's derivative of the order times the
directions' components there. The autoregressive coefficients are the
Levinson-Durbin map \\\xi\\ read at \\\rho = h(z)\\, with \\h\\ acting
componentwise. Faa di Bruno's formula over the set partitions \\\pi\\ of
the slots (the directions and the two free slots) gives \$\$D^n(\xi\circ
h) = \sum\_{\pi}\xi^{(\|\pi\|)}\big\[h^{(\|\beta\|)} \circ
\textstyle\prod\_{s\in\beta}s\big\]\_{\beta\in\pi},\$\$ since the inner
map is diagonal. A block holding no free slot pushes a vector, a block
holding one free slot scales that slot, and a block holding both puts a
diagonal. Each \\\xi^{(j)}\\ is read from
[`.gas_levinson_contract()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_levinson_contract.md)
against the pushed vectors.

## See also

[`gas_levinson5()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson5.md)
