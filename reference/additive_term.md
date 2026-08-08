# S7 Class for Additive Terms

The branch of
[`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
for terms that contribute a linear block \\X_j \beta_j\\ to the linear
predictor of one distribution parameter. A built additive term carries
the design block, the coefficient names, the blueprint that reproduces
the mapping on new data, and the penalty attached to its coefficients
(`NULL` when the term is unpenalized).

## Usage

``` r
additive_term(
  label = character(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL
)
```

## Arguments

- label:

  A character string prefixed to the coefficient names when non-empty.

- X:

  The design block, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- coef_names:

  The coefficient names, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- blueprint:

  The information needed to reproduce the mapping on new data, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- penalty:

  The penalty on the term's coefficients, or `NULL`.

## Value

An object inheriting from class `additive_term`.

## See also

[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)

## Examples

``` r
S7::S7_inherits(linpar(~1), additive_term)
#> [1] TRUE
```
