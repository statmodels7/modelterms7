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

Without the response every level starts at zero and every gap at one,
which puts the regimes wherever the scale of the data does not; on
response in the hundreds that start is far from any mode. The quantiles
spread the regimes over the data, as
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
spreads break-points over a covariate. `target` is the response on the
scale of the predictor, which the fitting layer supplies where the
family reads the parameter directly; where it is absent, or where two
quantiles coincide, the conventional start is returned.
