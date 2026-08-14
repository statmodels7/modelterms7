# S7 Class for Break-Point Terms

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a covariate whose effect changes at estimated break-points: a change
of slope
([`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md)), a
change of level
([`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)),
or both at the same points
([`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)).
The design block is the working one of the iteration that estimates the
break-points, and is recomputed by
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
as they move.

## Usage

``` r
SegTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  kind = character(0),
  var = NULL,
  npsi = integer(0),
  linear = logical(0),
  subformulas = list(),
  spec = list()
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

- kind:

  Which of the three constructions.

- var:

  The covariate expression.

- npsi:

  The number of break-points.

- linear:

  Whether the block carries the linear effect.

- subformulas:

  A named list of one-sided formulas, one per developed parameter.

- spec:

  The resolved construction settings.

## Value

An object of class `SegTerm`.

## See also

[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)

## Examples

``` r
S7::S7_inherits(seg(x), SegTerm)
#> [1] TRUE
```
