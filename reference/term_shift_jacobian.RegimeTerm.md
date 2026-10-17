# The Posterior Shift of a Regime Term, and Its Jacobian

Computes the shift and its Jacobian by forward mode through the
normalized forward and backward recursions of
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md),
over each group in time order.

## Arguments

- term:

  A built
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta, y, logdens, grad, psi, seed, cols, level:

  As for
  [`term_shift_jacobian()`](https://statmodels7.github.io/modelterms7/reference/term_shift_jacobian.md).

- ...:

  Unused.

## Value

A list with `shift` and `jacobian`.
