# Log-Likelihood of a Marginal Break-Point Term

Computes the exact marginal likelihood of a break-point model with
latent positions per group, together with the exact derivatives in the
term's own parameters on the parameter scale.

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
\\2^K\\ side patterns whose transition factors over the coordinates, and
each flip is weighted by its interval's prior mass, so the cost is \\n K
2^K\\, against the \\(n+1)^K\\ of the cells. The derivatives are
propagated through the same recursion: the derivatives of the interval
masses with respect to the parameters of the prior, and the derivatives
of the emissions with respect to the changes of level. For the
continuous kinds the conditional is smooth within an interval and the
integral runs on Gauss-Kronrod panels
([`numericals7::gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.html)).
The interior nodes lie at fixed points determined by the data, so the
derivatives in the parameters of the prior involve only the prior. The
region below the data, where the hinge keeps moving, is covered by
panels that follow the bulk of the prior, and the derivatives carry the
motion of their nodes. The region above the data, where the conditional
is constant, contributes its closed tail mass.

Every panel spans at most 2.5 prior standard deviations. The region
below the data carries at least four panels, and each interval between
two consecutive distinct observations at least one, so the number of
panels grows as `psi1.sd` decreases against the spacing of the
covariate. A group that would need more than 4000 panels is rejected.

The per-observation contributions are the one-step predictive densities
given the group's observations at smaller covariate values, and they sum
to each group's marginal log-likelihood. They are the increments of the
running logarithm of the forward pass that gives the total.
