# The Penalty Entries a Random-Effect Term Declares

One entry over the whole block when the prior reads a group's effects
together, and one per within-group column when they are independent.

## Usage

``` r
.random_entries(prior, d, m, wnames)
```

## Arguments

- prior:

  The effects' distribution, or a list of one per column.

- d:

  The number of within-group columns.

- m:

  The number of groups.

- wnames:

  The within-group column names.

## Value

A list of entries, as
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
documents.

## Details

The coefficients are ordered group by group, so column \\j\\ is the
stride \\j, d+j, 2d+j, \dots\\ – a subset of the term's own parameters,
NAMED rather than selected with a map, which is what keeps a kinked
prior's proximal operator available: a separable penalty under a
selection map is the generalized-lasso problem and has none.
