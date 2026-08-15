# S7 Class for Penalized Parametric Terms

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a parametric block whose coefficients carry a penalties7 penalty.
Constructed by
[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and [`mcp`](https://statmodels7.github.io/modelterms7/reference/mcp.md);
the four differ only in the penalty their factory attaches at build
time, and every derivative, hyperparameter, bound, link and kink is the
penalty object's.

## Usage

``` r
PenalizedTerm(
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
  input = NULL,
  input_expr = NULL,
  factory = function() NULL,
  standardize = FALSE
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

- input:

  The block as given: a one-sided formula or a numeric matrix.

- input_expr:

  The expression that produced a matrix input, kept so
  [`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  can re-evaluate it in new data.

- factory:

  The function mapping a coefficient count to the penalty object. It is
  called with the diagonal map as a second argument where `standardize`
  asks for one, so a factory that will never be standardized may take
  the count alone.

- standardize:

  Whether the block's columns are put on a common scale by the penalty's
  diagonal map.

## Value

An object of class `PenalizedTerm`.

## See also

[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md)

## Examples

``` r
S7::S7_inherits(ridge(~x), PenalizedTerm)
#> [1] TRUE
```
