# The Columns of a Component That Belong to Each of Its Sub-Terms

Splits a developed parameter's columns among the sub-terms developing
it, in the order the block binds them. It is what fills the `sub_index`
field of a
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
entry.

## Usage

``` r
component_sub_index(index, subs)
```

## Arguments

- index:

  An integer vector: the component's columns in the term's block, as
  long as the sub-terms' coefficient counts sum to.

- subs:

  A list of built sub-terms developing the parameter, in the order their
  blocks were bound.

## Value

A list of integer vectors, one per sub-term, partitioning `index` in
order. An empty list when `subs` is empty.

## Details

A developed parameter's block is its sub-terms' blocks bound side by
side in the order they were given, so the division is their coefficient
counts cumulated: with counts \\k_1, \dots, k_m\\ the \\i\\-th sub-term
takes `index` at positions \\k_1 + \dots + k\_{i-1} + 1\\ to \\k_1 +
\dots + k_i\\. The counts come from
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
so every sub-term must be built.

The term computes this rather than a consumer, because it rests on how
the block was assembled.

## See also

[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md),
the only caller;
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
for the counts it divides by.
