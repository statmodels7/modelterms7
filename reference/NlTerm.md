# S7 Class for Nonlinear Parametric Terms

The subclass of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
holding a parametric function that is nonlinear in its own parameters.
Its design block is the **Jacobian** of that function, so the term is
linear in the sense a fitting layer needs while the function is not. The
block depends on where the parameters currently are, and
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
recomputes it as they move.

## Usage

``` r
NlTerm(
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
  fn = NULL,
  nl_params = character(0),
  links = list(),
  subformulas = list(),
  deriv_mode = character(0),
  spec = list()
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

- ids:

  Which of the term's hyperparameters are shared with those of other
  terms, and under what label: a character vector named by the term's
  own hyperparameters, or `character(0)` for none. The terms carrying
  the same label for the same hyperparameter estimate one value. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md).

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

- fn:

  The function or formula defining the contribution, as given.

- nl_params:

  A character vector of the parameter names, in the order the derivative
  components are keyed by.

- links:

  A named list of one linkfunctions7 link per parameter.

- subformulas:

  A named list of one-sided formulas, one per parameter developed over
  covariates. Empty where none is.

- deriv_mode:

  `"symbolic"` or `"numeric"`.

- spec:

  A named list of the resolved construction settings.

## Value

An S7 object of class `NlTerm`, inheriting from
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the six properties above beside the ten they supply.

## The six properties of its own

`fn` is the function or formula as given, and `spec$is_formula` records
which. `nl_params` names the parameters, inferred from a formula or
taken from `params` for an opaque function.

`links` holds one linkfunctions7 link per parameter, the identity where
none was named, and `subformulas` one right-hand side per parameter
developed over covariates, empty otherwise.

`deriv_mode` is `"symbolic"` where
[`stats::deriv()`](https://rdrr.io/r/stats/deriv.html) could read the
expression and `"numeric"` where the derivatives are differenced.
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)
chooses per order, so the field records the route the expression as a
whole took.

`spec` carries the resolved construction settings, and the blueprint
what the build computed: the parameter designs, the starting
coefficients and the coefficients last committed.

## What makes it an additive term

While \\f\\ is differentiable the term is an ordinary additive one at
every point, so it needs no branch of its own in a fitting layer. What
it adds is three generics:
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
to recompute the block,
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
to report \\f\\ itself, and
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
to say that the block really is a Jacobian, which licenses a
Gauss-Newton step and a line search on the model's own objective.

## See also

[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md), the
constructor;
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)
for the derivatives;
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
and
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
for the pair a Gauss-Newton step reads.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0, 3, length.out = 60))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(60, sd = 0.05)

tm <- nl(~ a * exp(-r * x), start = list(a = 1, r = 1))
S7::S7_inherits(tm, NlTerm)
#> [1] TRUE
tm@nl_params
#> character(0)
tm@deriv_mode
#> [1] "symbolic"

# An opaque function cannot be read symbolically.
nl(function(x, theta) theta$a * exp(-theta$r * x),
   params = c("a", "r"), x = x, start = list(a = 1, r = 1))@deriv_mode
#> [1] "numeric"

# The block is the Jacobian, so a linear fit on it is a Gauss-Newton step.
b <- term_build(tm, dd)
dim(term_matrix(b))
#> [1] 60  2
term_jacobian_block(b)
#> [1] TRUE
```
