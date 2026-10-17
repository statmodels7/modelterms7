# Posterior Break-Points of a Marginal Term

Returns the posterior mean and standard deviation of each group's
break-points. For the step kind under the gaussian prior these are
mixtures of truncated-normal moments over the interval posterior, and
the moments of the edge intervals are computed through
[`numericals7::mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.html);
under an explicit prior the truncated moments come from
[`distributions7::truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.html)
and
[`distributions7::expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.html),
one interval at a time. For the continuous kinds they are the moments of
the node posterior, and the closed upper tail enters through its
truncated-normal moments.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

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
