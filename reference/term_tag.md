# The Covariance Label of a Term

The label a term carries to say that its coefficients share a covariance
block with those of other terms, or `NA_character_` where it carries
none.

## Usage

``` r
term_tag(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A single string, or `NA_character_` where the term carries no label.

## Details

The label is written in the middle position of a bar formula,
`random(~ 1 + x | u | id)`, following brms: the last position is the
grouping variable and the middle one the label. Terms carrying the same
label and the same grouping describe effects that are correlated with
each other, and a fitting layer collects them into one block.

The label is a **name and not data**. It is read as a symbol and never
evaluated, so a column of the data with the same name is not looked at.

The base method returns `NA_character_`, which is the answer for every
term but
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
so a term written outside this package needs no method here.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
the one constructor that reads a label;
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the penalties a term declares on its own.

## Examples

``` r
# No label: the effects of this term correlate with nothing outside it.
term_tag(random(~ 1 | g))
#> [1] NA
term_tag(linpar(~ x))
#> [1] NA

# A label, read from the middle of the two bars.
term_tag(random(~ 1 + x | u | g))
#> [1] "u"

# It is a name, not data: a column called `u` is not read.
d <- data.frame(x = rnorm(6), u = 1:6, g = factor(rep(c("a", "b"), 3)))
identical(term_coef_names(term_build(random(~ 1 | u | g), d)),
          term_coef_names(term_build(random(~ 1 | g), d)))
#> [1] TRUE
```
