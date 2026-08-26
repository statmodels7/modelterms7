# The Chart's Derivatives in the Unconstrained Parameters

The level, the score loadings and the autoregressive coefficients as
functions of the term's unconstrained parameters, with their first and
second derivatives in those.

## Usage

``` r
.gas_chart_derivs(zeta, p, q, links)
```

## Arguments

- zeta:

  The term's parameters on the unconstrained scale.

- p, q:

  The score and autoregressive orders.

- links:

  The links, as
  [`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
  gives them.

## Value

A list with the values and the derivative arrays.

## Details

The level and the loadings each reach the recursion through their own
link, so their first derivative is the link's and their second, on the
diagonal, is
[`linkfunctions7::d2linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d2linkinv.html);
on the identity both collapse to one and zero. The persistence reaches
the coefficients through two maps, the link onto the partial
autocorrelations and Levinson-Durbin onto the coefficients, so its
second derivative carries both a term in the map's own curvature and one
in the link's.
