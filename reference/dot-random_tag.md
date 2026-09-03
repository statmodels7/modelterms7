# The Label of a Random-Effect Term, Normalized

The term's `tag` as a single string, `NA_character_` where it carries
none.

## Usage

``` r
.random_tag(term)
```

## Arguments

- term:

  A
  [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md).

## Value

A single string, possibly `NA_character_`.

## Details

The property is a character vector and is empty on a term built before
the label existed, so the two spellings of "no label" are read in one
place rather than at each caller.

## See also

[`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md),
which returns it.
