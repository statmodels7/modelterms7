# The Score-Driven Recursion's Fourth Derivative

The fourth derivative of the filtered predictor in the caller's
unknowns, contracted against the caller's weights and against two
directions, and the lower orders' contractions beside it.

## Arguments

- term:

  A built `GasTerm`.

- eta:

  The static part of the predictor.

- y:

  The response, unused directly.

- score, curvature:

  The callbacks of
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The parameters on the parameter scale.

- g:

  The weights the fourth derivative is contracted against.

- seed:

  The derivative of the static predictor in the unknowns.

- blocks:

  The model's derivative pieces; see
  [`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md).

- directions:

  A list of two directions to contract against.

- ...:

  Unused.

## Value

A list with `jacobian`, `dphi`, `dpsi` and `curvature`.

## Details

The recursion carries five states – \\F_t\\, \\\Phi_t\\, the two third
derivatives \\\Psi^{(v)}\_t\\ and \\\Psi^{(w)}\_t\\, and \\\Omega_t\\ –
and each of the recursion's products is differentiated by
[`.gas_prod4()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_prod4.md),
the sixteen-term rule written once. The persistence's own fourth
derivative comes from
[`gas_levinson4()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson4.md)
through
[`.gas_chart_derivs4()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs4.md),
and the starting level's is read off the fixed point \\f_0 = \omega +
Sf_0\\.

Both third derivatives are needed and not one: the product rule at this
order pairs each direction's third derivative with the other direction.
