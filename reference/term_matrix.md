# Design Block of a Built Term

The \\n \times k\\ design block of a built additive term, with the
term's coefficient names as column names.

## Usage

``` r
term_matrix(term, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- ...:

  Passed to methods.

## Value

A numeric matrix.

## See also

[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md),
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`term_coef_names`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md),
[`term_npar`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
[`term_is_built`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)

## Examples

``` r
term_matrix(term_build(linpar(~x), data.frame(x = 1:4)))
#>   (Intercept) x
#> 1           1 1
#> 2           1 2
#> 3           1 3
#> 4           1 4
```
