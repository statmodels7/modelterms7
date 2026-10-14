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

A list with `within`, the left side of the bar (a call, a symbol or the
constant `1`), `group`, a symbol or a call, and `tag`, a single string
or `NA_character_`.

## Details

R nests the bars to the left, so in `~ 1 + x | u | id` the outer bar has
`id` on its right and `1 + x | u` on its left. The last position is
therefore the grouping variable and the middle one the label, which is
brms's convention and is the reading that the parser itself yields.

The label is taken with
[`as.character()`](https://rdrr.io/r/base/character.html) on the symbol
and **never evaluated**, so a column of the data carrying that name is
not read. It is the one element of the formula that does not refer to
the data.

Two bars are the most that carry a meaning, so a third is rejected.

The constructor and
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
both call this function, so that the split is defined in one place.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
and
[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md),
its two callers;
[`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md)
for what reports the label afterwards.
