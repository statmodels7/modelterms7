# S7 Class for Additive Terms

The branch of
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
whose terms contribute a block of design columns \\X_j \beta_j\\ to the
linear predictor of one distribution parameter. Everything a formula can
write that is not a filter or a latent mixture is on this branch:
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
the five penalized terms,
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md),
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
the break-point terms.

A built term on this branch carries four things beyond the root's
properties: the block itself, the coefficient names, the blueprint that
reproduces the block on new rows, and the penalty on its coefficients.

## Usage

``` r
additive_term(
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
  penalty = NULL
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

## Value

Nothing: the class is abstract and cannot be instantiated.
`additive_term()` throws. As a type it is the parent of every term that
contributes design columns, carrying the four properties above beside
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)'s
six.

## What a build fills, and what stays empty

A specification has `X` empty and `blueprint` an empty list;
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills all four. The blueprint is the reason
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
can **reapply** the mapping instead of rebuilding it: the factor levels,
the contrasts, the knots a basis was placed on and the spreads a
standardization used are recorded there at build time.
[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)'s
subset check is what tests that they really are.

`penalty` is `NULL` for an unpenalized term. Reading it directly answers
for a term carrying one penalty over its whole block;
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
is the general form and covers a term whose penalty reaches part of its
parameters or which carries several.

## The contract on this branch

[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
returns `X`,
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
returns `coef_names`,
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
counts the columns and
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
returns `penalty`, all from the properties, so a subclass gets them
without writing anything. What a subclass owes is
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
and
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md).

A term whose block moves with its own coefficients, as
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
do, also implements
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
and
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md),
and answers
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
to say whether the block it returns is a Jacobian or a frozen working
linearization.

## The block need not be dense

`X` is `class_any`, so a block may be a base matrix or any Matrix class.
A grouping indicator has one non-zero per row, and
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
builds it sparse;
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
and the penalized terms build sparse when asked. Code reading a block
tests `is.matrix(x) && is.numeric(x)` or a two-dimensional S4 object,
since [`is.matrix()`](https://rdrr.io/r/base/matrix.html) alone is
`FALSE` for every Matrix class.

## See also

[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
for the other branch,
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
and
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
for the contract, and
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
for the simplest term on this branch.

## Examples

``` r
d <- data.frame(x = rnorm(20), g = factor(rep(letters[1:4], 5)))

# A specification carries no block; a build fills all four properties.
spec <- ridge(~ x + g)
c(X = length(spec@X), blueprint = length(spec@blueprint))
#>         X blueprint 
#>         0         0 
b <- term_build(spec, d)
dim(term_matrix(b))
#> [1] 20  5
term_coef_names(b)
#> [1] "ridge.x"  "ridge.ga" "ridge.gb" "ridge.gc" "ridge.gd"
b@penalty
#> quadratic penalty on 5 coefficient(s) through 5 row(s); theta: lambda

# An unpenalized term has a NULL penalty.
is.null(term_build(linpar(~ x), d)@penalty)
#> [1] TRUE

# A grouping indicator is built sparse, so the block is not a base matrix.
class(term_matrix(term_build(random(~ 1 | g), d)))
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"

# Abstract: there is nothing to construct.
try(additive_term())
#> Error in S7::new_object(S7::S7_object(), label = label, hyper = hyper,  : 
#>   Can't construct an object from abstract class <additive_term>
```
