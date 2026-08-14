# Where a Nonlinear Term's Coefficients Begin

The coefficients
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
built the block at: the starting value of each of the term's own
parameters carried through its link, from `start` where the caller gave
one. The block is the Jacobian at those values, so it is not the same
block at any other point, and a fitting layer that started at zero would
linearize where the term was never meant to be evaluated.

## Arguments

- term:

  A built
  [`NlTerm`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).

- ...:

  Unused.

## Value

A numeric vector, one value per column of the block.
