# The Penalties of a Smooth Term

One entry over the whole block where a single roughness matrix covers
it, which is the default; one per level of a factor `by` under
`by_hyper = "level"`; and one over the penalized coordinates alone where
the smoother carries a penalty factory.

## Arguments

- term:

  A built
  [SmoothTerm](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md).

- ...:

  Unused.

## Value

A list of entries, as
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
documents.

## Details

The default answer is the base reading of
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
returned by calling that method rather than by repeating it, so the two
cannot drift.
