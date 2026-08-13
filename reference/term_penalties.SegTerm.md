# Penalties of a Segmented Term

One entry per kind of change the term carries a penalty on, naming the
coefficients it covers: `"delta"` for the slope changes, `"kappa"` for
the jump sizes, shared across the levels of `by`. The list is empty when
`penalty = NULL`, and for a specification, whose parameters there is
nothing yet to index: a penalty is attached at build, as it is for every
penalized term here.

## Arguments

- term:

  A built
  [`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- ...:

  Unused.

## Value

A list of entries, as
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
documents.
