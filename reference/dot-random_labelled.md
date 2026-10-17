# A Labelled Random-Effect Term, Built Without a Penalty

Finishes
[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md)
for a term carrying a covariance label: the block and its names as
usual, and no penalty entry at all.

## Usage

``` r
.random_labelled(term, Z, cn, parts, g, tt, mf, contr, d, wn)
```

## Arguments

- term:

  The term being built.

- Z:

  The block.

- cn:

  Its column names.

- parts:

  The formula's three pieces.

- g:

  The grouping factor.

- tt, mf, contr:

  The within-group terms, model frame and contrasts.

- d:

  The number of within-group columns.

- wn:

  Their names.

## Value

The built term, with no penalty.

## Details

The coefficients of a labelled term share a covariance block with those
of every other term carrying the same label and the same grouping, and
that block may span another equation entirely. The prior over it is
therefore a property of the class and not of any of its members, and a
penalty built here would give two priors for the same block, with no
rule for which the fit should read.
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
returns an empty list, and the fitting layer collects the labelled terms
and declares one entry over their stacked columns.

The blueprint records what the layer needs for that: the within-group
dimension and the grouping's levels, beside the pieces that
[`term_predict.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_predict.RandomTerm.md)
already reads.

The function rejects `hyper`, `id` and `correlated = FALSE` on a
labelled term, and a `distrib` that is a list or is univariate. A
`distrib` given here is the class's joint prior. Its dimension is the
total over the class and cannot be checked from one term, so that check
is left to the layer. The remaining checks, which verify that it is a
continuous distributions7 object whose location is held, are made here,
where the message can name the argument.

## See also

[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md),
its caller;
[`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md).
