# The Within-Group Design of a Term With Grouped Effects

The design row each group's effect multiplies, at new rows.

## Usage

``` r
term_within(term, newdata, ...)
```

## Arguments

- term:

  A built term.

- newdata:

  A data frame.

- ...:

  Passed to methods.

## Value

A numeric matrix of `nrow(newdata)` rows, or `NULL`.

## Details

A term whose coefficients are one effect per group,
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
builds its block by interacting a within-group design with the group
indicators. This returns the first factor alone, which does not depend
on the grouping: a caller integrating the effect of a group the fit
never saw reads \\z_i^\top b\\ from it with \\b\\ drawn from the prior.

The base method returns `NULL`, which is the answer for every term but
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md).

## See also

[`term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.md)
for the layout of the block,
[`term_within.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_within.RandomTerm.md)
for the method.

## Examples

``` r
d <- data.frame(x = rnorm(6), g = factor(rep(c("a", "b", "c"), 2)))
term_within(term_build(random(~ 1 + x | g), d), d[1:2, ])
#>   (Intercept)          x
#> 1           1 -1.2313234
#> 2           1  0.9838956
term_within(term_build(linpar(~ x), d), d)
#> NULL
```
