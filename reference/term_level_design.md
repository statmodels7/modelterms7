# The Design of a Term's Level Development

Where the parameter
[`term_level_param`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)
names is developed with covariates, the design of that development, with
one column per coordinate named as
[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
names it. `NULL` for a scalar level and for every other term.

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
`term_level_param` answers for the constant: a coordinate whose column
is constant shifts the equation's predictor exactly as an intercept
does. With the level developed, a direction of the development's span
that also lies in the span of the equation's design raises the same
question for that direction, and only a fitting layer, which holds both
designs, can ask it. This generic hands it the one half it cannot see.

## See also

[`term_level_param`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)

## Examples

``` r
is.null(term_level_design(linpar(~x)))
#> [1] TRUE
```
