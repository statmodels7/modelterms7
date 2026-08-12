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
  [`term_links`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
  gives them.

## Value

A list with the values and the derivative arrays.

## Details

The level and the loadings carry the identity link, so their first
derivative is one and their second is zero. The persistence reaches the
coefficients through two maps – the link onto the partial
autocorrelations and Levinson-Durbin onto the coefficients – so its
second derivative carries both a term in the map's own curvature and one
in the link's, which is where
[`d2linkinv`](https://statmodels7.github.io/linkfunctions7/reference/d2linkinv.html)
enters.
