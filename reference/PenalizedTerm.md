# S7 Class for Penalized Parametric Terms

The subclass of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a parametric block whose coefficients carry a penalties7 penalty.
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
all construct it; the five differ only in the penalty their factory
attaches at build time, so every derivative, hyperparameter, bound, link
and kink belongs to the penalty object and none of it is restated here.

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
  standardize = FALSE,
  sparse = NULL
)
```

## Arguments

- label:

  A character string prefixed to the term's coefficient names when
  non-empty, and used as the title of
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and the tag
  of [`print()`](https://rdrr.io/r/base/print.html). `character(0)` and
  `""` both mean no label.

- hyper:

  The hyperparameters of the term's penalty the caller **held**, as a
  named list keyed by the penalty's names. Empty, the default, means
  every one of them is estimated. See
  [`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md).

- grid:

  How many values a path visits for each of the term's hyperparameters,
  as a named list of single whole numbers. Empty, the default, leaves
  the number to the fitting layer. See
  [`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md).

- values:

  The values a path visits, for each hyperparameter the caller wrote
  out, as a named list of numeric vectors. Empty, the default, has the
  path build them. See
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

- X:

  The design block, one row per observation: a numeric matrix or a
  two-dimensional Matrix. Empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- coef_names:

  The block's coefficient names, one per column, prefixed by `label`
  when there is one. Filled by
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- blueprint:

  A named list of everything needed to reproduce the mapping on new
  rows. Its contents are the subclass's business; nothing outside the
  term reads them. Empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- penalty:

  A penalties7 penalty on the block's coefficients, or `NULL` when the
  term is unpenalized.

- input:

  The block as given: a one-sided formula or a numeric matrix.

- input_expr:

  The expression that produced a matrix input, kept so that
  [`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  can re-evaluate it in new data. `NULL` on the formula route.

- factory:

  A function of the coefficient count, returning the penalties7 penalty
  over that many coefficients, and taking the diagonal map as a second
  argument where `standardize` asks for one.

- standardize:

  A single logical: whether the block's columns are put on a common
  scale by the penalty's diagonal map. `FALSE` by default.

- sparse:

  `TRUE`, `FALSE` or `NULL` for the formula route's storage; see
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  for the rule `NULL` is settled by. A matrix input needs no such
  argument.

## Value

An S7 object of class `PenalizedTerm`, inheriting from
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the five properties above beside the ten they supply.

## The five properties of its own

`input` is the block as it was given, a one-sided formula or a matrix.
`input_expr` is the expression that produced a matrix input, kept so
that
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
can re-evaluate it in new data.

`factory` maps a coefficient count to the penalty. It is called at
build, when the count is finally known, and with the diagonal map as a
second argument when `standardize` asks for one, so a factory that will
never be standardized may take the count alone.

`standardize` is a single logical on the specification; after a build
the spreads themselves are in `blueprint$standardize`, one per column.

`sparse` governs the formula route only, a matrix input being kept in
whatever storage it arrives in. `NULL` until the build settles it.

## What the class does not decide

Which hyperparameters are estimated, over what grid and how far down a
path are recorded in
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)'s
own properties, read through
[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md),
[`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md),
[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md),
[`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md)
and
[`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md).
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
is `TRUE` for a built
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
and `FALSE` for the other four, read from each penalty's kink set.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five constructors share,
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
and its four siblings for their own hyperparameters,
[`print.PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/print.PenalizedTerm.md)
for how one displays.

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))

# Every one of the five returns this class.
vapply(list(ridge(~ x1), lasso(~ x1), enet(~ x1), scad(~ x1), mcp(~ x1)),
       function(t) S7::S7_inherits(t, PenalizedTerm), logical(1))
#> [1] TRUE TRUE TRUE TRUE TRUE

# The penalty is attached at build, when the column count is known.
spec <- lasso(~ x1 + x2)
c(spec = is.null(spec@penalty))
#> spec 
#> TRUE 
b <- term_build(spec, dd)
b@penalty
#> separable [fixed laplace2 [mu=0]] penalty on 2 coefficient(s) through 2 row(s); theta: lambda
b@penalty@n_coef
#> [1] 2

# standardize is a flag before the build and the spreads after it.
dd$x3 <- 1000 * dd$x2
bs <- term_build(lasso(~ x1 + x3, standardize = TRUE), dd)
bs@standardize
#> [1] TRUE
bs@blueprint$standardize
#>    lasso.x1    lasso.x3 
#>   0.7822996 886.9946868 
```
