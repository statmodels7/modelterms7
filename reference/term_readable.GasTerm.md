# What a Fitted Score-Driven Term Reports

The level, the score loadings and the autoregressive coefficients of the
literature, `omega`, `kappa1` and `xi1`, with the Jacobian from the
term's own parameters.

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

The level and the loadings are reported through their own links, each a
function of its own coordinate alone, which the base method already
does. The persistence is not: it is carried on a partial
autocorrelation, and the coefficients come from the Levinson-Durbin
recursion, whose Jacobian the term already computes for the filter.
Chained onto the rhobit link of each coordinate, that Jacobian is what a
delta-method standard error for \\\xi_j\\ needs. At \\q = 1\\ the two
coincide and the chain factor is the link's alone; above it they do not.

The scale an interval is built on follows the same split. At \\q = 1\\
it is the link of the partial autocorrelation (the rhobit unless `links`
gives another), so the interval stays inside the coordinate's range.
Above it the coefficients range over the stationary region, which is not
a box, and the interval is built on the identity scale.

A deviation is reported as it stands, being unconstrained and defined on
the scale of the parameter it departs from.
