# The Grouping by Which a Term's Coefficients Are Indexed

For a term whose coefficients are one set per level of a grouping
variable, returns the expression that defines that variable, its levels,
and the number of columns that each level carries. The result is `NULL`
for every other term.

## Usage

``` r
term_group(term, ...)
```

## Arguments

- term:

  A built term. An unbuilt one returns `NULL`.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A list with `expr` (a language object), `levels` (a character vector),
`dim` (a single integer) and `names` (`dim` column names), or `NULL`.

## Details

Two labeled terms belong to the same covariance block only if they share
a grouping, and deciding that requires both the expression and the
levels: the expression by itself is not enough, since `droplevels(id)`
and `id` are different expressions for the same grouping and `id` under
two subsets is the same expression for two different ones. The levels
settle it, and the expression makes a message readable.

`dim` is the number of columns that one level carries, which together
with the levels gives the layout of the block: the coefficients are
group-major, so level \\i\\ occupies positions \\(i-1)d + 1\\ to \\id\\,
where \\d\\ is `dim`. A consumer stacking two such blocks into one
covariance needs this layout, and `names` carries the within-group
design's own column names in the same order, so that each coordinate of
the covariance that it builds is identified as belonging to `x` or to
the intercept and is not merely numbered.

The base method returns `NULL`, which is the answer for every term but
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md).

## See also

[`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md)
for the label that, with this, identifies a covariance block;
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
for the formula from which both are read.

## Examples

``` r
d <- data.frame(x = rnorm(6), g = factor(rep(c("a", "b", "c"), 2)))
b <- term_build(random(~ 1 + x | g), d)
g <- term_group(b)
c(dim = g$dim, levels = length(g$levels))
#>    dim levels 
#>      2      3 
g$names
#> [1] "(Intercept)" "x"          

# every other kind of term has none
term_group(term_build(linpar(~ x), d))
#> NULL
```
