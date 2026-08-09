# Number of Coefficients of a Built Term

The number of columns of the term's design block.

## Usage

``` r
term_npar(term, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- ...:

  Passed to methods.

## Value

An integer.

## See also

[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md),
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`term_matrix`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md),
[`term_coef_names`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md),
[`term_is_built`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)

## Examples

``` r
term_npar(term_build(linpar(~x), data.frame(x = 1:4)))
#> [1] 2
```
