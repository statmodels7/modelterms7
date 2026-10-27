# S7 Class for the Unpenalized Parametric Term

`LinparTerm` is the subclass of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
that holds an unpenalized parametric block built from a one-sided
formula through
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
constructs it, and
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
collects a formula's bare covariates into it, so most models carry one
whether or not the caller wrote it.

Beyond the additive branch's four properties it records the formula,
whether the block is sparse, the contrasts used for the formula's
factors, and whether the column `(Intercept)` is dropped after the
coding.

## Usage

``` r
LinparTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  ids = character(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  formula = NULL,
  sparse = NULL,
  contrasts = list(),
  drop_intercept = FALSE
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

- blueprint:

  A named list of everything needed to reproduce the mapping on new
  rows. Its contents belong to the subclass, and only the term reads
  them. It is empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- penalty:

  A penalties7 penalty on the block's coefficients, or `NULL` when the
  term is unpenalized.

- formula:

  The one-sided formula defining the block, such as `~ x + log(z) + f`.
  Its environment is kept and used when a symbol is absent from the
  data.

- sparse:

  `TRUE` to build the block as a `dgCMatrix` through
  [`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
  `FALSE` for a base matrix, `NULL` to let the build choose. It is kept
  as given; the storage actually used is in `blueprint$sparse`.

- contrasts:

  A named list of contrasts for the formula's factors, in
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)'s
  own form, or an empty list for the session's defaults.

- drop_intercept:

  `TRUE` to remove the column `(Intercept)` after the model matrix is
  built, so that the factors keep the coding that they have in a formula
  with an intercept.
  [`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
  sets it on a
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  written in a formula that already carries an intercept. The default is
  `FALSE`.

## Value

An S7 object of class `LinparTerm`, inheriting from
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the four properties above beside the properties that it inherits.

## The properties of its own

`formula` is the one-sided formula, kept with its environment, so a
symbol that the data do not carry is still found where the formula was
written.

`sparse` holds the value that the caller gave, and `NULL` where the
caller gave no value. The build leaves it unchanged and records in
`blueprint$sparse` the storage that it resolved, which is the value that
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reads, so the block and a prediction from it never differ in storage. A
formula naming a factor of many levels has one non-zero per row, and
building the dense model matrix of it is the memory use that sparse
storage avoids.

`contrasts` is a named list, one entry per factor, or empty for the
session's own `options("contrasts")`. Whatever is used is recorded in
the blueprint and reapplied, so a fit and a prediction always use the
same coding.

`drop_intercept` is `FALSE` by default. When it is `TRUE`, the build
removes the column `(Intercept)` after the coding, and the blueprint
records the setting so that a prediction applies it again.

The class carries no penalty: `penalty` is `NULL` on every `LinparTerm`,
and
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
returns an empty list, so
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
counts its columns exactly.

## See also

[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
the constructor to use;
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
which builds one implicitly;
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
and
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md).

## Examples

``` r
d <- data.frame(x = rnorm(20), g = factor(rep(letters[1:4], 5)))

# linpar() is the constructor and returns an object of this class.
tm <- linpar(~ x + g)
S7::S7_inherits(tm, LinparTerm)
#> [1] TRUE
tm@formula
#> ~x + g
#> <environment: 0x5631d3e1ea18>

# The property keeps the value given; the blueprint records the storage
# that the build used, and a prediction reads the blueprint.
bt <- term_build(tm, d)
c(sparse_is_null = is.null(bt@sparse), settled = bt@blueprint$sparse)
#> sparse_is_null        settled 
#>           TRUE          FALSE 
term_build(linpar(~ g, sparse = TRUE), d)@blueprint$sparse
#> [1] TRUE

# The contrasts are recorded and reapplied, so the coding cannot drift.
b <- term_build(linpar(~ g, contrasts = list(g = "contr.sum")), d)
term_coef_names(b)
#> [1] "(Intercept)" "g1"          "g2"          "g3"         

# It is never penalized.
c(penalty = is.null(b@penalty), entries = length(term_penalties(b)),
  edf = edf(b))
#> penalty entries     edf 
#>       1       0       4 
```
