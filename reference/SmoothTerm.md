# S7 Class for Smooth Terms

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a penalized smooth of one or more covariates: a basis7 expansion
with a roughness penalty on its coefficients. Constructed by
[`s`](https://statmodels7.github.io/modelterms7/reference/s.md) for one
covariate and by
[`te`](https://statmodels7.github.io/modelterms7/reference/te.md) for
several.

## Usage

``` r
SmoothTerm(
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
  vars = list(),
  by = NULL,
  spec = list(),
  sparse = FALSE
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

- vars:

  The expressions of the covariates being smoothed.

- by:

  An optional expression the smooth varies with.

- spec:

  The construction settings: the basis, its dimension and degree, and
  whether the linear part is carried separately.

- sparse:

  Whether the block is a `dgCMatrix`, which only a factor `by` admits.
  See [`s`](https://statmodels7.github.io/modelterms7/reference/s.md).

## Value

An object of class `SmoothTerm`.

## See also

[`s`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te`](https://statmodels7.github.io/modelterms7/reference/te.md)

## Examples

``` r
S7::S7_inherits(s(x), SmoothTerm)
#> [1] TRUE
```
