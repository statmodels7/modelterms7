# Where a Break-Point Term's Coefficients Begin

The start
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
computed: the break-points at the positions `psi` names or at the
interior quantiles of the covariate, and unit changes. Where `target` is
given (the response on the scale of the predictor) and no coefficient
carries a development, the slope and the changes are instead the
least-squares coefficients of `target` on an intercept and the term's
exact columns at the starting positions, which puts the starting mean
near the data whatever the scale of the covariate. Zero is degenerate
here, never neutral: a discontinuous term reads its break-point off
\\-g_k/\delta_k\\, which at zero is the same clamped position for every
one of them, and a continuous term's Jacobian column vanishes, so a
fitting layer that starts every coefficient at zero has to be told
otherwise.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- target:

  The response on the scale of the predictor, one value per observation
  of the build data, or `NULL` for the unit changes.

- ...:

  Unused.

## Value

A numeric vector, one value per column of the block.
