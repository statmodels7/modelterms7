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
  hyper = list(),
  grid = list(),
  min_ratio = numeric(0),
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

- hyper:

  The hyperparameters of the term's penalty that the caller HELD, as a
  named list. Empty, the default, means every one of them is estimated.
  See
  [`term_hyper`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md).

- grid:

  How many values a path visits for each of the term's hyperparameters,
  as a named list. Empty, the default, leaves it to the criterion. See
  [`term_grid`](https://statmodels7.github.io/modelterms7/reference/term_grid.md).

- min_ratio:

  How far down the path over the size of the kink reaches, as a fraction
  of the value that empties the block, or `numeric(0)` for the
  criterion's own. See
  [`term_path_min`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md).

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
