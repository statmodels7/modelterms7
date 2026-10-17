# Penalties of a Break-Point Term

Returns one entry per penalty that the developments of the term's own
coefficients declare, and each entry names the coefficients that it
covers. The entries of a subformula shared by every coefficient of one
kind are pooled into one, so that `gamma ~ 0 + lasso(~1)` is a single
penalized block over the changes under one hyperparameter. The list is
empty for a term whose coefficients carry no development, and for a
specification, which has no coefficients to index yet: the penalties are
recorded at build, as for every penalized term in the package.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- ...:

  Unused.

## Value

A list of entries, as
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
documents.
