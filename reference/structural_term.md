# S7 Class for Structural Terms

The branch of
[`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
for terms that rewrite the likelihood contributions rather than adding a
design block:
[`gas`](https://statmodels7.github.io/modelterms7/reference/gas.md),
whose predictor is a recursion, and
[`regime`](https://statmodels7.github.io/modelterms7/reference/regime.md),
whose contribution is a likelihood mixed over latent states. A term on
this branch reports its own parameters through
[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
rather than coefficients, so
[`term_npar`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
counts those and
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
indexes into them.

## Usage

``` r
structural_term(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0)
)
```

## Arguments

- label:

  A character string naming the term.

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

## Value

An object inheriting from class `structural_term`.

## See also

[`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md)

## Examples

``` r
S7::S7_inherits(linpar(~1), structural_term)
#> [1] FALSE
```
