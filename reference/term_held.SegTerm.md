# Coefficients Held by a Break-Point Term

Returns the coefficient slots that carry the break-points of a term in
the held state set by
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md).
For a
[`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
that is not held, built or not, the result is an empty vector.

## Arguments

- term:

  A
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- ...:

  Unused.

## Value

An integer vector of positions in the term's coefficients.
