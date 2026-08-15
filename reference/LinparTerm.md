# S7 Class for the Unpenalized Parametric Term

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for the unpenalized parametric block built from a one-sided formula
through [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html).
Constructed by
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
and implicitly by
[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
which collects the bare covariates of a model formula into one term of
this class.

## Usage

``` r
LinparTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  formula = NULL
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

- search:

  How the term's own hyperparameters are covered when it has several
  with a kink: `"grid"` for every combination of them, `"cyclic"` for
  one at a time, or `character(0)` for the default. See
  [`term_search`](https://statmodels7.github.io/modelterms7/reference/term_search.md).

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

- formula:

  The one-sided formula defining the block.

## Value

An object of class `LinparTerm`.

## See also

[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md)

## Examples

``` r
S7::S7_inherits(linpar(~1), LinparTerm)
#> [1] TRUE
```
