# The Covariance Label of a Random-Effect Term

The label written between two bars, or `NA_character_` where the formula
carries one bar.

## Arguments

- term:

  A
  [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md),
  built or not.

- ...:

  Unused.

## Value

A single string, or `NA_character_`.

## Details

It is recorded by the constructor, from the middle position of
`~ lhs | tag | g`, and is a name that is never evaluated.

## See also

[`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md)
for the generic,
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
for the formula it is read from.
