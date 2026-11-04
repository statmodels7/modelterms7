# S7 Class for Shape-Constrained Tensor Product Smooths

`ScteTerm` is the subclass of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a smooth of two covariates with a shape (monotone, convex or
concave, or a combination) in one or both of them:
[`scte()`](https://statmodels7.github.io/modelterms7/reference/scte.md)
constructs it. The tensor-product coefficients are written through the
maps of the margins' shapes, so the contribution is nonlinear in the
free coefficients \\\gamma\\ and the block is its Jacobian at the
current coefficients, as for
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md).

## Usage

``` r
ScteTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  ids = character(0),
  X = NULL,
  coef_names = character(0),
  build_info = list(),
  penalty = NULL,
  vars = list(),
  by = NULL,
  shapes = list(),
  spec = list(),
  sparse = NULL
)
```

## Arguments

- label:

  A character string prefixed to the term's coefficient names when
  non-empty, and used as the title of
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and the tag
  of [`print()`](https://rdrr.io/r/base/print.html). The empty string
  `""` means no label.

- hyper:

  The hyperparameters of the term's penalty that the caller held, as a
  named list keyed by the penalty's names. Empty, the default, means
  every one of them is estimated. See
  [`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md).

- grid:

  How many values a path visits for each of the term's hyperparameters,
  as a named list of single whole numbers. Empty, the default, leaves
  the number to the fitting layer. See
  [`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md).

- values:

  The values that a path visits, for each hyperparameter that the caller
  wrote out, as a named list of numeric vectors. Empty, the default, has
  the path build them. See
  [`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md).

- min_ratio:

  How far down the path over the size of the kink reaches, as a fraction
  of the value that empties the block: one number in \\(0, 1)\\, or
  `numeric(0)` for the fitting layer's own. See
  [`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md).

- search:

  How the term's own hyperparameters are covered when it has more than
  one carrying a kink: `"grid"` for every combination of them,
  `"cyclic"` for one at a time, or `character(0)` for the default. See
  [`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md).

- ids:

  Which of the term's hyperparameters are shared with those of other
  terms, and under what label: a character vector named by the term's
  own hyperparameters, or `character(0)` for none. The terms carrying
  the same label for the same hyperparameter estimate one value. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md).

- X:

  The design block, one row per observation: a numeric matrix or a
  two-dimensional Matrix. It is empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- coef_names:

  The block's coefficient names, one per column, prefixed by `label`
  when there is one.
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- build_info:

  A named list of everything needed to reproduce the mapping on new
  rows. Its contents belong to the subclass, and only the term reads
  them. It is empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- penalty:

  A penalties7 penalty on the block's coefficients, or `NULL` when the
  term is unpenalized.

- vars:

  A list of the two covariate expressions, unevaluated.

- by:

  The `by` expression, unevaluated, or `NULL`.

- shapes:

  A list of two elements, each a monotone or curvature
  [shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md)
  or `NULL` for a free margin.

- spec:

  A named list of the settings that the build uses: the `smoothers`,
  `default_label` and `by_hyper`.

- sparse:

  `FALSE`.

## Value

An S7 object of class `ScteTerm`, inheriting from
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md).

## See also

[`scte()`](https://statmodels7.github.io/modelterms7/reference/scte.md),
the constructor.

## Examples

``` r
S7::S7_inherits(scte(x, z), ScteTerm)
#> [1] TRUE
```
