# The Columns of a Component That Belong to Each of Its Sub-Terms

Splits a developed parameter's columns among the sub-terms developing
it, in the order in which the block binds them. The result fills the
`sub_index` field of a
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
entry.

## Usage

``` r
component_sub_index(index, subs)
```

## Arguments

- index:

  An integer vector: the component's columns in the term's block, whose
  length is the sum of the sub-terms' coefficient counts.

- subs:

  A list of built sub-terms developing the parameter, in the order in
  which their blocks were bound.

## Value

A list of integer vectors, one per sub-term, partitioning `index` in
order. The result is an empty list when `subs` is empty.

## Details

A developed parameter's block is its sub-terms' blocks bound side by
side in the order in which they were given, so the division is their
coefficient counts cumulated: with counts \\k_1, \dots, k_m\\ the
\\i\\-th sub-term takes `index` at positions \\k_1 + \dots + k\_{i-1} +
1\\ to \\k_1 + \dots + k_i\\. The counts come from
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
so every sub-term must be built.

The term computes this division, and a consumer does not, because it
rests on how the block was assembled.

## See also

[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md),
whose methods are its only callers;
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
for the counts by which it divides the columns.
