# Log-Likelihood of a Marginal Break-Point Term

The exact marginal likelihood of a break-point model with latent
positions per group, with the exact derivatives in the term's own
parameters on the parameter scale.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- eta:

  The static predictor.

- y:

  The response, reaching the sum through the callbacks.

- logdens, score:

  The log-density and its derivative in the predictor.

- psi:

  The parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Unused.

## Value

A list with `loglik` and `jacobian`, the latter on the parameter scale.

## Details

For the step kind the conditional is constant on the product partition
of the intervals between a group's ordered observations, and the exact
sum over it is taken by the forward recursion of the side chain: the
monotone process of active break-points is a hidden Markov chain on the
\\2^K\\ side patterns whose transition factors over the coordinates,
each flip weighted by its interval's prior mass, so the cost is \\n K
2^K\\, against the \\(n+1)^K\\ of the cells. The derivatives ride the
same recursion: the masses' in the prior's parameters, the emissions' in
the changes of level. For the continuous kinds the conditional is smooth
within an interval and the integral runs on a fixed Gauss-Kronrod panel
per interval
([`numericals7::gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.html)),
the interior nodes fixed points of the data so that the derivatives in
the prior's parameters read the prior alone; the region below the data,
where the hinge keeps moving, is covered by panels that follow the
prior's bulk, whose node motion the derivatives carry, and the region
above it, where the conditional is constant, contributes its closed tail
mass.

The per-observation contributions are the one-step predictive densities
given the group's observations at smaller covariate values, which sum to
each group's marginal log-likelihood; producing them costs one pass over
the component weights per observation, so the decomposition is a factor
of the group's size dearer than the total alone, which the per-group
cell or node count already prices.
