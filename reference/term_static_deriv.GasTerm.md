# The Filter's Sensitivity to Its Own Equation's Coefficients

The forward recursion of
[`term_static_deriv`](https://statmodels7.github.io/modelterms7/reference/term_static_deriv.md)
for a score-driven term.

## Arguments

- term:

  A built score-driven term.

- curv:

  The curvature at each predictor.

- X:

  The directions to propagate.

- psi:

  The term's parameters, on the parameter scale.

- ...:

  Ignored.

## Value

A matrix of `X`'s dimensions.

## Details

The recursion is the filter's own with the seed replaced: the starting
level depends on the term's parameters and not on the static predictor,
so the propagated derivative starts at zero and everything after it
comes from the scores. The autoregressive and loading coefficients are
read at each row, which covers a term whose parameters carry submodels
of their own as well as one whose parameters are constant.

It is written in R rather than compiled: it evaluates no callback, runs
once at a variance rather than at every iteration of a fit, and costs
one pass over the data per column.

## See also

[`term_static_deriv`](https://statmodels7.github.io/modelterms7/reference/term_static_deriv.md),
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
