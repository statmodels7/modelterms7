# The Design of a Term's Level Development

Returns, where the parameter that
[`term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)
names is developed with covariates, the design of that development, with
one column per coordinate named as
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
names it. The value is `NULL` for a scalar level and for every other
term.

## Usage

``` r
term_level_design(term, ...)
```

## Arguments

- term:

  A built term.

- ...:

  Passed to methods.

## Value

A numeric matrix with named columns, or `NULL`.

## Details

It exists for the subspace form of the confounding question.
`term_level_param` covers the case of a constant: a coordinate whose
column is constant shifts the equation's predictor exactly as an
intercept does. If the level is developed, a direction of the
development's span that also lies in the span of the equation's design
raises the same question for that direction. The question can be settled
only in a fitting layer, which holds both designs, and this generic
supplies the design of the development, which the layer cannot compute
from the term alone.

## See also

[`term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)
for the scalar case,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for the names that its columns carry.

## Examples

``` r
is.null(term_level_design(linpar(~x)))
#> [1] TRUE
```
