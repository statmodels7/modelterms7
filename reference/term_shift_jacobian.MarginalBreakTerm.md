# The Posterior Shift of a Marginal Break-Point Term, and Its Jacobian

Computes the shift and its Jacobian. For the step kind, the function
applies forward mode through the forward and backward recursions of the
side chain over the \\2^K\\ side patterns, and the interval masses and
the survivals carry their derivatives in the prior's parameters. For the
continuous kinds, it computes the derivative of each group's posterior
over its quadrature nodes, and the node weights, the node positions and
the shifts carry their derivatives in the term's own parameters.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- eta, y, logdens, grad, psi, seed, cols, level:

  As for
  [`term_shift_jacobian()`](https://statmodels7.github.io/modelterms7/reference/term_shift_jacobian.md).

- ...:

  Unused.

## Value

A list with `shift` and `jacobian`.
