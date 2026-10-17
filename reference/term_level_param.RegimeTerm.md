# The Level of a Regime Term

Returns `"level1"`. The levels are ordered by construction: each is the
previous one plus a positive gap, so the first one shifts every regime
and is the direction that an intercept in the same equation also spans.
A constant leaves the gaps unchanged, because it cannot express a
difference between regimes.

## Arguments

- term:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- ...:

  Unused.

## Value

A single string.
