# The Chart's Fourth Derivatives, in Two Directions

[`.gas_chart_derivs3()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs3.md)
differentiated once more and contracted against a second direction: one
matrix for the level, one per score loading and one per autoregressive
coefficient.

## Usage

``` r
.gas_chart_derivs4(zeta, p, q, links, vz, wz)
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

- vz, wz:

  The two directions, in the same coordinates as `zeta`.

## Value

A list with `q_omega`, `q_a` and `q_b`, each a matrix or a list of
matrices over the term's base coordinates.

## Details

The level and the loadings ride scalar links, so each of their fourth
derivatives is a single diagonal entry, the link's own \\h''''\\ times
the two directions' components there. The persistence is a composition,
the Levinson-Durbin map read at \\\rho = h^{-1}(z)\\, whose inner map is
DIAGONAL: differentiating \\B(z) = \phi(\rho(z))\\ four times and
contracting two slots leaves, with \\p_v = h'v\\ and \\p_w = h'w\\ the
directions pushed onto the partial autocorrelations,
\$\$Q\[p_v,p_w\]h'h'^{\top} + T\[p_v\]\big(h''w\\h'^{\top} +
h'\\h''w^{\top}\big) + T\[p_w\]\big(h''v\\h'^{\top} +
h'\\h''v^{\top}\big) + T\[h''vw\]h'h'^{\top}\$\$ \$\${}+
H\big(h'''vw\\h'^{\top} + h'\\h'''vw^{\top} + h''v\\h''w^{\top} +
h''w\\h''v^{\top}\big) + \mathrm{diag}\big(h''T\[p_v\]p_w +
h''H(h''vw) + h'''w\\Hp_v + h'''v\\Hp_w + h''''vw\\P\big),\$\$ with
\\Q\\, \\T\\, \\H\\ and \\P\\ the map's fourth, third, second and first
derivatives. The expression is symmetric in \\v\\ and \\w\\, which a
caller can read off it and a test asserts.

## See also

[`.gas_chart_derivs3()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs3.md),
[`gas_levinson4()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson4.md)
