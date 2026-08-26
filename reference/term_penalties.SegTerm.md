# Penalties of a Break-Point Term

One entry per penalty the developments of the term's own coefficients
declare, naming the coefficients it covers. The entries of a subformula
shared by every coefficient of one kind are pooled into one, so that
`gamma ~ 0 + lasso(~1)` is a single penalized block over the changes
under one hyperparameter. The list is empty for a term whose
coefficients carry no development, and for a specification, whose
coefficients there is nothing yet to index: a penalty is attached at
build, as it is for every penalized term here.

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
