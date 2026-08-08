# S7 Class for the Unpenalized Parametric Term

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for the unpenalized parametric block built from a one-sided formula
through [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html).
Constructed by
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
and implicitly by
[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
which collects the bare covariates of a model formula into one term of
this class.

## Usage

``` r
LinparTerm(
  label = character(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  formula = NULL
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

- formula:

  The one-sided formula defining the block.

## Value

An object of class `LinparTerm`.

## See also

[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md)

## Examples

``` r
S7::S7_inherits(linpar(~1), LinparTerm)
#> [1] TRUE
```
