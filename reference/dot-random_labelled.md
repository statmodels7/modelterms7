# A Labelled Random-Effect Term, Built Without a Penalty

Finishes
[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md)
for a term carrying a covariance label: the block and its names as
usual, and no penalty entry at all.

## Usage

``` r
.random_labelled(term, Z, cn, parts, g, tt, mf, contr, d)
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

## Value

The built term, with no penalty.

## Details

The coefficients of a labelled term share a covariance block with those
of every other term carrying the same label and the same grouping, and
that block may span another equation entirely. The prior over it is
therefore a property of the class and not of any of its members, and
building one here would leave two answers to the same question with
nothing to say which the fit should read.
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
returns an empty list, which is the right answer rather than a gap, and
the fitting layer collects the labelled terms and declares one entry
over their stacked columns.

What the blueprint records is what the layer needs to do that: the
within-group dimension and the grouping's levels, beside the pieces
[`term_predict.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_predict.RandomTerm.md)
already reads.

A `distrib` given here names the **class's** joint prior. Its dimension
is the total over the class and cannot be checked from one term, so that
one check is left to the layer; the rest – that it is a continuous
distributions7 object with its location held – is asked here, where the
message can name the argument.

## See also

[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md),
its caller;
[`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md).
