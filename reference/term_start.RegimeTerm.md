# Where a Regime Term's Parameters Start

The levels at the quantiles of the response at \\(2j-1)/(2k)\\, \\j = 1,
\dots, k\\, and the transition matrix with every row uniform.

## Arguments

- term:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- ...:

  Unused.

- target:

  The response on the scale of the predictor, or `NULL`.

## Value

A named numeric vector on the unconstrained scale.

## Details

Without the response the start is the zero vector on the unconstrained
scale: the first level is zero, every gap is one and every row of the
transition matrix is uniform, so the levels are \\0, 1, \dots, k - 1\\
whatever the scale of the data. The quantiles spread the regimes over
the data, as
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
spreads break-points over a covariate. `target` is the response on the
scale of the predictor, which the fitting layer supplies where the
family reads the parameter directly. Where it is absent, where fewer
than \\k\\ of its values are finite, or where two quantiles coincide,
the conventional start above is returned.
