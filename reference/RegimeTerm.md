# S7 Class for Markov Regime Terms

A subclass of
[`structural_term`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
for a latent Markov chain of regimes, each shifting the linear predictor
by a level of its own. Constructed by
[`regime`](https://statmodels7.github.io/modelterms7/reference/regime.md).

## Usage

``` r
RegimeTerm(
  label = character(0),
  k = integer(0),
  by = NULL,
  time = NULL,
  chain = NULL,
  blueprint = list()
)
```

## Arguments

- label:

  A character string prefixed to the term's coefficient names when
  non-empty.

- k:

  The number of regimes.

- by:

  An optional grouping expression, run independently.

- time:

  An optional ordering expression.

- chain:

  The parameters7 transition matrix.

- blueprint:

  The resolved ordering and grouping.

## Value

An object of class `RegimeTerm`.

## See also

[`regime`](https://statmodels7.github.io/modelterms7/reference/regime.md)

## Examples

``` r
S7::S7_inherits(regime(2), RegimeTerm)
#> [1] TRUE
```
