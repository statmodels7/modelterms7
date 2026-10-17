# S7 Class for Nonlinear Parametric Terms

`NlTerm` is the subclass of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
that holds a parametric function that is nonlinear in its own
parameters. Its design block is the Jacobian of that function, so the
term is linear in the sense that a fitting layer needs while the
function is not. The block depends on the current values of the
parameters, and
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

- fn:

  The function or formula defining the contribution, as given.

- nl_params:

  A character vector of the parameter names, in the order used to key
  the derivative components.

- links:

  A named list of linkfunctions7 links, one per parameter once the term
  is built.

- subformulas:

  A named list of one-sided formulas, one per parameter developed over
  covariates. The list is empty if the term develops no parameter over
  covariates.

- deriv_mode:

  One of `"symbolic"`, `"formula_fd"` or `"numeric"`.

- spec:

  A named list of the resolved construction settings.

## Value

An S7 object of class `NlTerm`, inheriting from
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the properties described above in addition to those of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md).

## Properties of its own

`fn` is the function or formula as given, and `spec$is_formula` records
which. `nl_params` names the parameters:
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) sets
it from `params` for a function, and
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
sets it for a formula, from the names that the formula uses and that the
data do not supply.

`links` holds the links given in the call. After
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
it holds one linkfunctions7 link per parameter, the identity link for
each parameter without a named link. `subformulas` holds one right-hand
side per parameter developed over covariates, and is empty otherwise.

`deriv_mode` is `"symbolic"` where
[`stats::deriv()`](https://rdrr.io/r/stats/deriv.html) could read the
expression of a formula, `"formula_fd"` where it could not and the first
derivatives of the formula are differenced, and `"numeric"` for a
function, whose derivatives are differenced. A specification built from
a formula carries `"symbolic"` until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
sets the value. The property records the route of the first derivatives;
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)
chooses the route for each order separately.

`spec` carries the resolved construction settings. The blueprint carries
what the build computed: the parameter designs, the current coefficients
(the starting ones until
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
stores others), the value of the function at them, and the orders of
derivative that are analytic.

## What makes it an additive term

While \\f\\ is differentiable the term is an ordinary additive one at
every point, so a fitting layer needs no branch of its own for it. The
class implements
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
which recomputes the block,
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md),
which reports \\f\\ itself, and the derivatives of the block in the
coefficients
([`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md),
[`term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.md)
and
[`term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.md)).
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
is inherited from
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
and returns `TRUE`, which states that the block is a Jacobian and
licenses a Gauss-Newton step and a line search on the model's own
objective.

## See also

[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md), the
constructor;
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)
for the derivatives;
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
and
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
for the pair that a Gauss-Newton step reads.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0, 3, length.out = 60))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(60, sd = 0.05)

tm <- nl(~ a * exp(-r * x), start = list(a = 1, r = 1))
S7::S7_inherits(tm, NlTerm)
#> [1] TRUE

# The parameters are read from the formula when the term is built.
b <- term_build(tm, dd)
b@nl_params
#> [1] "a" "r"
b@deriv_mode
#> [1] "symbolic"

# An opaque function cannot be read symbolically.
nl(function(x, theta) theta$a * exp(-theta$r * x),
   params = c("a", "r"), x = x, start = list(a = 1, r = 1))@deriv_mode
#> [1] "numeric"

# The block is the Jacobian, so a linear fit on it is a Gauss-Newton step.
dim(term_matrix(b))
#> [1] 60  2
term_jacobian_block(b)
#> [1] TRUE
```
