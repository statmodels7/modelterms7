# The Grouping a Term's Coefficients Are Indexed By

For a term whose coefficients are one set per level of a grouping
variable, the expression that variable came from, the levels it took,
and how many columns each level carries. `NULL` for every other term.

## Usage

``` r
term_group(term, ...)
```

## Arguments

- term:

  A built term. An unbuilt one returns `NULL`.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A list with `expr` (a language object), `levels` (a character vector)
and `dim` (a single integer), or `NULL`.

## Details

Two labelled terms belong to the same covariance block only if they
share a grouping, and deciding that needs both halves: the expression
alone is not enough, since `droplevels(id)` and `id` are different
expressions for the same grouping and `id` under two subsets is the same
expression for two different ones. The levels settle it, and the
expression makes a message readable.

`dim` is the number of columns one level carries, which with the levels
says how the block is laid out: the coefficients are group-major, so
level \\i\\ occupies positions \\(i-1)d + 1\\ to \\id\\. A consumer
stacking two such blocks into one covariance needs exactly that.

The base method returns `NULL`, which is the answer for every term but
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md).

## See also

[`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md)
for the label that, with this, identifies a covariance block;
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
for the formula both are read from.

## Examples

``` r
d <- data.frame(x = rnorm(6), g = factor(rep(c("a", "b", "c"), 2)))
b <- term_build(random(~ 1 + x | g), d)
g <- term_group(b)
c(dim = g$dim, levels = length(g$levels))
#>    dim levels 
#>      2      3 

# every other kind of term has none
term_group(term_build(linpar(~ x), d))
#> NULL
```
