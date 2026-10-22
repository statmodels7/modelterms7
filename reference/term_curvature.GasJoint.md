# The Curvature of a Joint Score-Driven Filter

Computes the Jacobian of every filtered predictor in a caller's unknowns
and the second derivative of the filtered predictors contracted against
the caller's weights, propagated through the joint recursion.

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- eta, y, score, curvature:

  Unused: the recursion reads the data in `blocks`. They are kept for
  the signature of the generic.

- psi:

  The terms' parameters, as for
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- g:

  The weights, a matrix with one column per filtered parameter.

- seed:

  A list with one matrix per distribution parameter: the derivative of
  that equation's static predictor in the unknowns, one row per
  observation. The joint filter's parameters are the last columns, in
  the order of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- blocks:

  A list with `filtered` (the position in `seed` of the equation of each
  filtered parameter), `u` (the driving scores at the fitted predictors,
  a matrix with one column per filtered parameter), `U1` and `U2` (the
  first and second derivatives of each driving score in every predictor,
  arrays of dimension `c(n, K, P)` and `c(n, K, P, P)`), and, for the
  third and fourth orders, `U3` and `U4`, with one more dimension of
  length `P` each.

- ...:

  Unused.

## Value

A list with `jacobian` (one matrix per filtered parameter) and
`curvature`.

## Details

The level of each parameter and each driving score are carried as jets:
a value with its derivatives. The jet of a level is a sum of products of
the chart's jets with the jets of earlier driving scores and levels, and
the jet of a driving score at a row is formed from the jets of the
predictors of every equation by the chain rule, with the driving score's
derivatives in those predictors supplied as data in `blocks`. The static
rows of every equation enter through `seed`, and only a filtered
equation's predictor has a second derivative.
