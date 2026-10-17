# Drawing a Regime Path

Draws a path of the latent chain and returns the predictor of each
observation, which is the static predictor plus the level of the state
that the observation occupies.

## Arguments

- term:

  A built regime term.

- psi:

  The term's parameters, on the parameter scale.

- eta:

  The static part of the predictor, one value per row of the data on
  which the term was built.

- draw:

  Ignored: the levels do not depend on the response.

- ...:

  Ignored.

## Value

A list with `eta`, the static predictor plus the level of each
observation's state; `y` (`NULL`); and `latent`, the state of each
observation, in the order of the data.

## Details

One chain is drawn for each group of the blueprint, in time order, and
each starts from the chain's stationary distribution.
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
and
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
read the data in this way, and the likelihood is written with that
initial law, so starting from a fixed first state or from a uniform draw
would simulate a different model from the one that is fitted. The
transition matrix comes from the term's own
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html),
and the stationary law is the solution of \\\delta(I - P) = 0\\ with
\\\sum_j \delta_j = 1\\, which the method solves itself.

The response is not drawn, because the levels do not depend on it, so
`y` is returned as `NULL` and the caller draws at the returned
predictor. An `eta` whose length differs from the number of rows of the
build data signals an error.

## See also

[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md),
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
