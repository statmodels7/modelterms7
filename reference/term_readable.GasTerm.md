# What a Fitted Score-Driven Term Reports

Returns the level, the score loadings and the autoregressive
coefficients of the literature, `omega`, `kappa1` and `xi1`, with the
Jacobian from the term's own parameters.

## Arguments

- term:

  A
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).

- zeta:

  The parameters on the unconstrained scale.

- ...:

  Unused.

## Value

A list, as
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
documents.

## Details

The level and the loadings are reported through their own links, and
each of them is a function of only its own coordinate, as in the base
method. The persistence is not: it is carried on partial
autocorrelations, and the coefficients come from the Levinson-Durbin
recursion, which the term already differentiates for the filter. Each
partial autocorrelation is mapped through its own link, the rhobit
unless `links` gives another. The Jacobian of the recursion is chained
onto these links, and a delta-method standard error for \\\xi_j\\ needs
that chain. At \\q = 1\\ the two coincide and the chain factor reduces
to the link's derivative; above it they differ.

The scale on which an interval is built follows the same split. At \\q =
1\\ it is the link of the partial autocorrelation, so the interval stays
inside the coordinate's range. Above it the coefficients range over the
stationary region, which is not a box, and the interval is built on the
identity scale.

The coefficients of a development are reported as they stand, because
they are unconstrained and act on the scale of the link of the developed
parameter. If any partial autocorrelation is developed, the
autoregressive coefficients vary by observation, so a single \\\xi_j\\
does not exist, and the persistence coordinates are reported as the base
method gives them.
