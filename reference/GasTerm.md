# S7 Class for Score-Driven Dynamics

A subclass of
[`structural_term`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
for a generalized autoregressive score component: a time-varying level
driven by the score of the observation density, added to the linear
predictor. Constructed by
[`gas`](https://statmodels7.github.io/modelterms7/reference/gas.md).

## Usage

``` r
GasTerm(
  label = character(0),
  p = integer(0),
  q = integer(0),
  by = NULL,
  time = NULL,
  deviations = NULL,
  penalty_kind = NULL,
  blueprint = list()
)
```

## Arguments

- label:

  A character string prefixed to the term's coefficient names when
  non-empty.

- p:

  The number of score lags.

- q:

  The number of autoregressive lags.

- by:

  An optional grouping expression, filtered independently.

- time:

  An optional ordering expression.

- deviations:

  Which parameters carry a deviation per group.

- penalty_kind:

  The penalty on the deviations, if any.

- blueprint:

  The resolved ordering and grouping.

## Value

An object of class `GasTerm`.

## See also

[`gas`](https://statmodels7.github.io/modelterms7/reference/gas.md)

## Examples

``` r
S7::S7_inherits(gas(), GasTerm)
#> [1] TRUE
```
