# S7 Class for Grouped Random-Effect Terms

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for grouped coefficients with a distribution on the effects: the
within-group design interacted with the grouping indicators, one
coefficient per group and per within-group column, with the penalty
carrying the effects' distribution. Constructed by
[`random`](https://statmodels7.github.io/modelterms7/reference/random.md).

## Usage

``` r
RandomTerm(
  label = character(0),
  hyper = list(),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  formula = NULL,
  correlated = logical(0),
  precision = NULL,
  distrib = NULL,
  kinks = integer(0)
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

  The bar formula, e.g. `~ 1 | g` or `~ x | g`.

- correlated:

  Logical; whether the default Gaussian lets the within-group effects
  correlate.

- precision:

  A parameters7 matrix parameter for the per-group precision, or `NULL`.

- distrib:

  A univariate distributions7 object for the effects, or `NULL`.

- kinks:

  The declared kink set of `distrib`'s log-density.

## Value

An object of class `RandomTerm`.

## See also

[`random`](https://statmodels7.github.io/modelterms7/reference/random.md)

## Examples

``` r
S7::S7_inherits(random(~ 1 | g), RandomTerm)
#> [1] TRUE
```
