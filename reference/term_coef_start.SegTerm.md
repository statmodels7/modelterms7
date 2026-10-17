# The Starting Coefficients of a Break-Point Term

Returns the start that
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
computed: the break-points at the positions that `psi` names or at the
interior quantiles of the covariate, and unit changes. Where `target` is
given (the response on the scale of the predictor) and the coefficients
carry no development, the slope and the changes are instead the
least-squares coefficients of `target` on an intercept and the term's
exact columns at the starting positions, which puts the starting mean
near the data whatever the scale of the covariate. A start of zero is
degenerate for this term. A sharp discontinuous term reads its
break-point off \\-g_k/\delta_k\\ (a sharp jseg off a quadratic that
reduces to it), which at zero is the same clamped position for every
break-point, and a continuous term's Jacobian column vanishes, so a
fitting layer starts the term from this vector and not from zeros.

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
