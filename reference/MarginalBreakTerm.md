# S7 Class for Marginal Break-Point Terms

A subclass of
[`structural_term`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
for break-points that vary by group as latent variables integrated out
of the likelihood. Constructed by
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md) or
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
with `marginal = TRUE`.

## Usage

``` r
MarginalBreakTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  kind = character(0),
  var = NULL,
  npsi = integer(0),
  linear = logical(0),
  group = NULL,
  prior = NULL,
  spec = list(),
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

- values:

  The values a path visits, for each hyperparameter the caller wrote
  out, as a named list. Empty, the default, has the path build them. See
  [`term_values`](https://statmodels7.github.io/modelterms7/reference/term_values.md).

- min_ratio:

  How far down the path over the size of the kink reaches, as a fraction
  of the value that empties the block, or `numeric(0)` for the
  criterion's own. See
  [`term_path_min`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md).

- search:

  How the term's own hyperparameters are covered when it has several
  with a kink: `"grid"` for every combination of them, `"cyclic"` for
  one at a time, or `character(0)` for the default. See
  [`term_search`](https://statmodels7.github.io/modelterms7/reference/term_search.md).

- kind:

  Which of the three constructions.

- var:

  The covariate expression.

- npsi:

  The number of break-points.

- linear:

  Whether the term carries the linear effect as its own parameter (`seg`
  and `jseg`).

- group:

  The grouping expression, from the break-point's
  [`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
  subformula.

- prior:

  The latent's distribution: `NULL` for the gaussian, or a
  distributions7 object from `random(distrib = )`.

- spec:

  The resolved construction settings.

- blueprint:

  The resolved grouping and interval structure.

## Value

An object of class `MarginalBreakTerm`.

## See also

[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)

## Examples

``` r
S7::S7_inherits(jump(x, psi ~ random(~1 | id), marginal = TRUE),
                MarginalBreakTerm)
#> [1] TRUE
```
