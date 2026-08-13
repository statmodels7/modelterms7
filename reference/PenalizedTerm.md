# S7 Class for Penalized Parametric Terms

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a parametric block whose coefficients carry a penalties7 penalty.
Constructed by
[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
and
[`mcp`](https://statmodels7.github.io/modelterms7/reference/ridge.md);
the four differ only in the penalty their factory attaches at build
time, and every derivative, hyperparameter, bound, link and kink is the
penalty object's.

## Usage

``` r
PenalizedTerm(
  label = character(0),
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
