# What a Fitted Score-Driven Term Reports

The level, the score loadings and the AUTOREGRESSIVE COEFFICIENTS of the
literature – `omega`, `alpha1`, `beta1` – with the Jacobian from the
term's own parameters.

## Arguments

- term:

  A
  [`GasTerm`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).

- zeta:

  The parameters on the unconstrained scale.

- ...:

  Unused.

## Value

A list, as
[`term_readable`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
documents.

## Details

The level and the loadings are reported through their own links, each a
function of its own coordinate alone, which the base method already
does. The persistence is not: it is carried on a partial
autocorrelation, and the coefficients come from the Levinson-Durbin
recursion, whose Jacobian the term already computes for the filter.
Chained onto the rhobit link of each coordinate, that Jacobian is what a
delta-method standard error for \\\beta_j\\ needs. At \\q = 1\\ the two
coincide and the chain factor is the link's alone; above it they do not.

A deviation is reported as it stands, being unconstrained and defined on
the scale of the parameter it departs from.
