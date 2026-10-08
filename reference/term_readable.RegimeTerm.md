# What a Fitted Regime Term Reports

The first level, the gaps between consecutive levels and the transition
probabilities \\p\_{ij} = P(S_t = j \mid S\_{t-1} = i)\\, with the
Jacobian from the term's own parameters.

## Arguments

- term:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- zeta:

  The parameters on the unconstrained scale.

- ...:

  Unused.

## Value

A list, as
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
documents.

## Details

The levels and the gaps are reported through their own links, as the
base method does. The chain is estimated on the additive log-ratios of
each row of the transition matrix, which are coordinates and not the
quantities a reader wants, so every entry of the matrix is reported
instead, \\k^2\\ of them, row by row and named `p1.1`, `p1.2`, and so
on. Entry \\p\_{ij}\\ depends on the log-ratios of row \\i\\ alone, and
its row of the Jacobian is the derivative of
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html)'s
value. Each interval is built on the logit scale and mapped back, so it
stays inside \\(0, 1)\\.
