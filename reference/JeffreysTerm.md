# S7 Class for the Jeffreys Term

`JeffreysTerm` is the subclass of
[`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
that
[`jfp()`](https://statmodels7.github.io/modelterms7/reference/jfp.md)
constructs. It adds the power of the determinant and the setting of the
constant direction; the input handling, the storage and the prediction
are those of the penalized terms.

## Usage

``` r
JeffreysTerm(
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
  input = NULL,
  input_expr = NULL,
  penalty_fun = function() NULL,
  standardize = FALSE,
  sparse = NULL,
  power = integer(0),
  intercept = NULL
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

- input:

  The block as given: a one-sided formula or a numeric matrix.

- input_expr:

  The expression that produced the input, kept so that
  [`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  can re-evaluate a matrix input in new data. It is recorded on the
  formula route as well, where it is not read.

- penalty_fun:

  A function of the coefficient count, returning the penalties7 penalty
  over that many coefficients, and taking the diagonal map as a second
  argument when `standardize = TRUE`.

- standardize:

  A single logical: whether each coefficient is penalized on the scale
  of its own column, through the penalty's diagonal map. The default is
  `FALSE`.

- sparse:

  `TRUE`, `FALSE` or `NULL` for the formula route's storage; see
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  for the rule by which `NULL` is settled. A matrix input needs no such
  argument.

- power:

  The power of the determinant.

- intercept:

  The setting of the constant direction.

## Value

An S7 object of class `JeffreysTerm`, inheriting from
[`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md).

## Details

The build rejects a block whose columns are linearly dependent, by the
rank of the pivoted QR decomposition with the tolerance `1e-7` with
which statmodels7 detects aliased coefficients: the determinant of the
block's information would be zero.

## See also

[`jfp()`](https://statmodels7.github.io/modelterms7/reference/jfp.md)

## Examples

``` r
t <- jfp(~ x)
S7::S7_inherits(t, JeffreysTerm)
#> [1] TRUE
S7::S7_inherits(t, PenalizedTerm)
#> [1] TRUE
```
