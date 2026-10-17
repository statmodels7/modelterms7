# Whether a Break-Point Term's Break-Points Have Settled

Applies
[`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
through the generic, so that a fitting layer can test the stopping rule
of the iteration without inspecting the class of the term. A term that
has not been refreshed has not moved and reports `FALSE`, because the
first refresh is the one that measures nothing. A term held by
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md)
reports `TRUE`, because its positions are fixed.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- ...:

  Unused.

## Value

A single logical.
