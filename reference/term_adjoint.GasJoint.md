# Filter Several Score-Driven Terms Backwards

Runs the joint recursion of
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
in reverse and returns the derivative of a caller's objective with
respect to the static predictors of the filtered parameters and with
respect to the driving scores.

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- eta, y, score, curvature, psi:

  As for
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- g:

  The direct derivative of the objective in the filtered predictors, a
  matrix like `eta`.

- ...:

  Unused.

- stored_curv:

  `NULL`, or the `curv` that
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
  returned at the same parameters and static predictors; where it is
  `NULL`, the forward pass is run to obtain it. The name is not a prefix
  of `curvature`, so that it cannot be matched partially to that
  argument.

- fast, threads:

  Passed to that forward pass.

## Value

A list with `deta` and `dscore`, matrices like `eta`.

## Details

With \\\bar f\_{k,t}\\ the derivative of the objective in the level and
\\\bar u\_{l,t}\\ in the driving score, the objective reaches the
filtered predictor \\e\_{k,t}\\ directly through `g` and through every
driving score read there, \\\bar e\_{k,t} = g\_{k,t} + \sum_l \bar
u\_{l,t} C\_{t,lk}\\; the level passes this on to its own earlier values
through the persistence and to the earlier driving scores through the
loadings. The pass reads the curvature \\C_t\\ that the forward pass
stored and runs compiled.
