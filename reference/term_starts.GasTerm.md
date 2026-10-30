# The Starting Points of a Score-Driven Term

Returns `n_start` starting points. The first is
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)'s;
the following ones put the first partial autocorrelation of the
persistence at \\0.9\\, \\0.5\\, \\0.97\\ and \\0.7\\, in this order,
the other parameters starting as in the first.

## Arguments

- term:

  A
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).

- ...:

  Unused.

- target:

  Unused.

## Value

A list of `term@n_start` named numeric vectors on the unconstrained
scale; one vector when `n_start` is empty, 1, or the term has no
persistence.
