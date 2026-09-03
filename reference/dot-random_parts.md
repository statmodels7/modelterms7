# The Three Pieces of a Bar Formula

Splits `~ lhs | g` or `~ lhs | tag | g` into the within-group design,
the grouping expression and the covariance label, and rejects the shapes
that are not one of the two.

## Usage

``` r
.random_parts(formula)
```

## Arguments

- formula:

  The bar formula as given.

## Value

A list with `within` and `group`, both language objects, and `tag`, a
single string or `NA_character_`.

## Details

R nests the bars to the left, so in `~ 1 + x | u | id` the outer bar has
`id` on its right and `1 + x | u` on its left. The last position is
therefore the grouping variable and the middle one the label, which is
brms's convention and comes out of the parser rather than being imposed.

The label is taken with
[`as.character()`](https://rdrr.io/r/base/character.html) on the symbol
and **never evaluated**, so a column of the data carrying that name is
not read. It is the one thing in a formula that does not refer to the
data, and reading it as data is what the two-bar form did before: with
no such column it stopped with `object 'u' not found`, with a factor
column it built twice the columns and a warning about `'|'` applied to
factors, and with a numeric column it built twice the columns and said
nothing at all.

Two bars are the most that carry a meaning, so a third is rejected
rather than read.

Written once because the constructor and
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
both need the split and two copies of it would agree only by accident.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
and
[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md),
its two callers;
[`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md)
for what reports the label afterwards.
