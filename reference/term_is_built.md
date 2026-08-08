# Whether a Term Has Been Built

`TRUE` for a term returned by
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
and `FALSE` for a bare specification. The accessors
[`term_matrix`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md),
[`term_npar`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
[`term_coef_names`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
and
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reject a specification, and this predicate is the test they use.

## Usage

``` r
term_is_built(term)
```

## Arguments

- term:

  An object inheriting from class
  [`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md).

## Value

A logical scalar.

## Examples

``` r
term_is_built(linpar(~x))
#> [1] FALSE
term_is_built(term_build(linpar(~x), data.frame(x = 1:4)))
#> [1] TRUE
```
