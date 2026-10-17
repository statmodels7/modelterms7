# The Filter's Sensitivity to Its Own Equation's Coefficients

Computes the forward recursion of
[`term_static_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_static_deriv.md)
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

The recursion is the filter's own, except that the seed is replaced: the
starting level depends on the term's parameters, never on the static
predictor, so the propagated derivative starts at zero and everything
after it comes from the scores. The autoregressive and loading
coefficients are read at each row, which covers a term whose parameters
carry submodels of their own as well as one whose parameters are
constant.

It is written in R, evaluates no callback, and makes one pass over the
data, carrying all the columns of `X` together.

## See also

[`term_static_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_static_deriv.md),
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
