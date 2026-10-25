# A Static Parameter Whose Score Drives a Joint Filter

The S7 class of a member of a
[`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md)
that carries no level: the predictor of its parameter is static, and its
score enters the driving vector so that the levels of the other members
can load it through their `p_cross`.
[`gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.md)
creates it from its argument `drivers`.

## Usage

``` r
GasDriver(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  ids = character(0),
  blueprint = list(),
  p = integer(0),
  q = integer(0),
  by = NULL,
  time = NULL,
  links = list(),
  submodels = list(),
  scaling = numeric(0),
  p_cross = integer(0),
  cross_with = character(0)
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

- blueprint:

  A named list of the resolved ordering, grouping and sub-term designs,
  empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- p:

  The number of score lags, an integer of at least 1.

- q:

  The number of autoregressive lags, an integer of at least 0.

- by:

  An optional grouping expression; each group is filtered independently.
  `NULL` means one series.

- time:

  An optional ordering expression. `NULL` means row order.

- links:

  A named list of linkfunctions7 links that override the defaults, empty
  when the caller overrides no default.

- submodels:

  A named list of one-sided formulas, one per parameter developed over
  covariates. The list is empty if no parameter is developed.

- scaling:

  The power \\d\\ of the expected information in the factor
  \\\mathcal{I}^{-d}\\ that multiplies the score; zero for the score as
  it is.

- p_cross:

  A named integer vector with the number of lags of the score of each
  other filtered parameter that the level of this term loads, or an
  empty vector for a term that loads its own score only.

- cross_with:

  The distribution parameters whose scores the level of this term loads
  besides its own, in the order of their loadings. It is empty until
  [`gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.md)
  fills it from `p_cross`.

## Value

An S7 object of class `GasDriver`, inheriting from
[`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).

## Details

A driver is a score-driven term with no parameters, whose level is zero
at every observation:
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md),
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md),
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
and
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
are empty, and the recursion gives it no level and no derivative.

## See also

[`gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = rnorm(40))
jd <- gas_joint(list(mu = term_build(gas(p = 1, q = 1, time = t,
                                         p_cross = c(sigma = 1)), dd)),
                drivers = "sigma")
S7::S7_inherits(jd@terms$sigma, GasDriver)
#> [1] TRUE
term_params(jd@terms$sigma)
#> character(0)
```
