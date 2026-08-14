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
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL
)
```

## Arguments

- label:

  A character string prefixed to the coefficient names when non-empty.

- hyper:

  The hyperparameters of the term's penalty that the caller HELD, as a
  named list. Empty, the default, means every one of them is estimated.
  See
  [`term_hyper`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md).

- grid:

  How many values a path visits for each of the term's hyperparameters,
  as a named list. Empty, the default, leaves it to the criterion. See
  [`term_grid`](https://statmodels7.github.io/modelterms7/reference/term_grid.md).

- values:

  The values a path visits, for each hyperparameter the caller wrote
  out, as a named list. Empty, the default, has the path build them. See
  [`term_values`](https://statmodels7.github.io/modelterms7/reference/term_values.md).

- min_ratio:

  How far down the path over the size of the kink reaches, as a fraction
  of the value that empties the block, or `numeric(0)` for the
  criterion's own. See
  [`term_path_min`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md).

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
