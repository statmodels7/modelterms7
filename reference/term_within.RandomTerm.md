# The Within-Group Design of a Random-Effect Term

The design one group's effect multiplies, at new rows: the intercept,
the covariates of a random slope, one column per coordinate of the
effect.

## Arguments

- term:

  A built
  [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md).

- newdata:

  A data frame carrying the within-group covariates.

- ...:

  Unused.

## Value

A numeric matrix of `nrow(newdata)` rows and `term_group(term)$dim`
columns, named as `term_group(term)$names`.

## Details

Read from the blueprint the build recorded, so a factor keeps its levels
and contrasts. It does not depend on the grouping, which is what lets a
caller integrate the effect of a group the fit never saw: the effect is
\\z_i^\top b\\ with \\z_i\\ this row and \\b\\ drawn from the term's
prior.

## See also

[`term_within()`](https://statmodels7.github.io/modelterms7/reference/term_within.md)
for the generic,
[`term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.md)
for the layout.

## Examples

``` r
dd <- data.frame(x = c(0.5, 1, 2), g = factor(c("a", "b", "a")))
b <- term_build(random(~ 1 + x | g), dd)
term_within(b, data.frame(x = c(3, 4), g = c("new", "a")))
#>   (Intercept) x
#> 1           1 3
#> 2           1 4
```
