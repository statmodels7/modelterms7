# The Chart's Third Derivatives, in One Direction

[`.gas_chart_derivs`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs.md)
differentiated once more and contracted against a single direction in
the term's own coordinates: one matrix for the level, one per score
loading and one per autoregressive coefficient.

## Usage

``` r
.gas_chart_derivs3(zeta, p, q, links, vz)
```

## Arguments

- zeta:

  The term's base parameters on the unconstrained scale.

- p, q:

  The score and autoregressive orders.

- links:

  The links, as
  [`term_links`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
  gives them.

- vz:

  The direction, in the same coordinates as `zeta`.

## Value

A list with `t_omega`, `t_a` and `t_b`, each a matrix or a list of
matrices over the term's base coordinates.

## Details

The level and the loadings ride scalar links, so each of their third
derivatives is a single diagonal entry, the link's own \\h'''\\ times
the direction's component there. The persistence is a composition, the
Levinson-Durbin map read at \\\rho = h^{-1}(z)\\, and differentiating
\\B(z) = \phi(\rho(z))\\ three times and contracting the last slot gives
\$\$T\_{kl}h'\_kh'\_l + P\_{kl}(h''\_kv_kh'\_l + h'\_kh''\_lv_l) +
\delta\_{kl}\big((Hw)\_k h''\_k + P_kh'''\_kv_k\big),\$\$ with \\w_m =
h'\_mv_m\\ the direction pushed onto the partial autocorrelations, \\T\\
the contracted third derivative of the map and \\P\\, \\H\\ its first
two.

## See also

[`.gas_chart_derivs`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs.md),
[`gas_levinson3`](https://statmodels7.github.io/modelterms7/reference/gas_levinson3.md)
