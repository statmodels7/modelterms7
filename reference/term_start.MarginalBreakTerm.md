# Where a Marginal Break-Point Term Starts

Computes the starting values from a target (the response on the scale of
the predictor, which the fitting layer supplies) by a least-squares
profile in two stages. In the first stage, common positions are
selected, with a free intercept for each group, and the linear effect
and the changes of level and slope are estimated on all groups pooled.
In the second stage, the positions of each group are moved one at a
time, while the pooled coefficients are held fixed. Without a target,
the starting positions are quantiles of the covariate, the scale is
derived from its spread, and the changes are zero.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- ...:

  Unused.

- target:

  The response on the predictor's scale, or `NULL`.

## Value

A named numeric vector on the unconstrained scale.

## Details

With a target, the candidate positions are 23 quantiles of the covariate
between its 6% and 94% points. The local minima of the pooled residual
sum of squares over the candidates (at most four) are the seeds. For \\K
\> 1\\ the further positions are added greedily to a seed. Each group's
positions are then moved, one at a time, over the midpoints between its
distinct observations, and the seed whose final fit leaves the smallest
residual sum of squares is retained. That fit has coefficients shared
across groups and the positions of each group, and it supplies the
linear effect and the changes.

The medians of the group positions give `psik.mean`, unless positions
were given to the constructor, which then take precedence. Their median
absolute deviations give the prior scale, which is floored in every case
at half the median gap between consecutive distinct covariate values.
Without a target the positions are the quantiles of the covariate at \\k
/ (K + 1)\\ and the scale is one sixth of the distance between its 5%
and 95% quantiles. Under an explicit prior, `psi1.mean` takes the
position, a parameter of the prior that is interpreted as a scale starts
at the scale through its own link, and the other parameters start at
zero. A target of the wrong length is not used, and an unbuilt term
returns zeros.
