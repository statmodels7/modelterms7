# Coefficient Names of a Built Term

The names of the term's coefficients, prefixed by the term's label when
the label is non-empty.

## Usage

``` r
term_coef_names(term, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- ...:

  Passed to methods.

## Value

A character vector.

## See also

[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md),
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`term_matrix`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md),
[`term_npar`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
[`term_is_built`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)

## Examples

``` r
term_coef_names(term_build(linpar(~x), data.frame(x = 1:4)))
#> [1] "(Intercept)" "x"          
```
