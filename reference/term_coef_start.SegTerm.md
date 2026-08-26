# Where a Break-Point Term's Coefficients Begin

The start
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
computed: unit changes, and the break-points at the positions `psi`
names or at the interior quantiles of the covariate. Zero is degenerate
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

  Unused: a break-point term already reads the covariate's interior
  quantiles at
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- ...:

  Unused.

## Value

A numeric vector, one value per column of the block.
