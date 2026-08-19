# Posterior Break-Points of a Marginal Term

The posterior mean and standard deviation of each group's break-points.
For the step kind under the gaussian prior these are mixtures of
truncated-normal moments over the interval posterior, the edge intervals
read through
[`mills_ratio`](https://statmodels7.github.io/numericals7/reference/mills_ratio.html);
under an explicit prior the truncated moments come from
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.html)
and
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.html),
one interval at a time. For the continuous kinds they are the moments of
the node posterior, the closed upper tail entering through its
truncated-normal moments.

## Arguments

- term:

  A built
  [`MarginalBreakTerm`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens:

  The log-density.

- psi:

  The term's parameters.

- ...:

  Unused.

## Value

A data frame with `group`, `psi`, `mean` and `sd`.
