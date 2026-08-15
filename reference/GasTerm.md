# S7 Class for Score-Driven Dynamics

A subclass of
[`structural_term`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
for a generalized autoregressive score component: a time-varying level
driven by the score of the observation density, added to the linear
predictor. Constructed by
[`gas`](https://statmodels7.github.io/modelterms7/reference/gas.md).

## Usage

``` r
GasTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  p = integer(0),
  q = integer(0),
  by = NULL,
  time = NULL,
  links = list(),
  submodels = list(),
  blueprint = list()
)
```

## Arguments

- label:

  A character string prefixed to the term's coefficient names when
  non-empty.

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

- search:

  How the term's own hyperparameters are covered when it has several
  with a kink: `"grid"` for every combination of them, `"cyclic"` for
  one at a time, or `character(0)` for the default. See
  [`term_search`](https://statmodels7.github.io/modelterms7/reference/term_search.md).

- p:

  The number of score lags.

- q:

  The number of autoregressive lags.

- by:

  An optional grouping expression, filtered independently.

- time:

  An optional ordering expression.

- links:

  The links overriding the defaults, if any.

- submodels:

  One optional subformula per parameter.

- blueprint:

  The resolved ordering and grouping.

## Value

An object of class `GasTerm`.

## See also

[`gas`](https://statmodels7.github.io/modelterms7/reference/gas.md)

## Examples

``` r
S7::S7_inherits(gas(), GasTerm)
#> [1] TRUE
```
