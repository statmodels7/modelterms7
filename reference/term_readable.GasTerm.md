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
function of its own coordinate alone, as the base method does. The
persistence is not: it is carried on partial autocorrelations, and the
coefficients come from the Levinson-Durbin recursion, which the term
already differentiates for the filter. Each partial autocorrelation is
mapped through its own link, the rhobit unless `links` gives another,
and the Jacobian of the recursion is chained onto these links, which is
what a delta-method standard error for \\\xi_j\\ needs. At \\q = 1\\ the
two coincide and the chain factor is the link's alone; above it they do
not.

The scale an interval is built on follows the same split. At \\q = 1\\
it is the link of the partial autocorrelation, so the interval stays
inside the coordinate's range. Above it the coefficients range over the
stationary region, which is not a box, and the interval is built on the
identity scale.

The coefficients of a development are reported as they stand, being
unconstrained and acting on the scale of the link of the developed
parameter. If any partial autocorrelation is developed, the
autoregressive coefficients vary by observation, no single \\\xi_j\\
exists, and the persistence coordinates are reported as the base method
gives them.
