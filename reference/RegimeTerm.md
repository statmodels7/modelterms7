# S7 Class for Markov Regime Terms

The subclass of
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
holding a latent Markov chain of regimes, each shifting the linear
predictor by a level of its own.
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
constructs it. Its contribution is a likelihood mixed over the
unobserved state path, so it implements
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
and has no predictor to report.

## Usage

``` r
RegimeTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  blueprint = list(),
  k = integer(0),
  by = NULL,
  time = NULL,
  chain = NULL
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

- blueprint:

  A named list of the resolved ordering and grouping, empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- k:

  The number of regimes, an integer of at least 2.

- by:

  An optional grouping expression; each group runs its own recursion.
  `NULL` for one group.

- time:

  An optional ordering expression. `NULL` for row order.

- chain:

  A
  [`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html)
  of side `k`.

## Value

An S7 object of class `RegimeTerm`, inheriting from
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the five properties above beside
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)'s
six.

## The five properties of its own

`k` is the number of regimes. `by` and `time` are the grouping and
ordering expressions as written, kept unevaluated; `NULL` means one
group in row order.

`chain` is the
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html)
object, whose free values are the additive log-ratios of each row, so
every row is a probability vector at any coordinate. Its `free_names`
are the tail of
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

`blueprint` is filled by
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
and holds the row order within each group and the observation count. The
class overrides the branch's `blueprint` property because a structural
term carries no design block to hang one on.

## See also

[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md),
the constructor;
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
for what it computes;
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
for the smoothed states;
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
for the other dynamic term.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))

tm <- regime(2, time = t)
S7::S7_inherits(tm, RegimeTerm)
#> [1] TRUE
c(k = tm@k, chain = class(tm@chain)[1])
#>                                    k                                chain 
#>                                  "2" "parameters7::TransitionMatrixParam" 

# The build resolves the order and fills the blueprint.
b <- term_build(tm, dd)
names(b@blueprint)
#> [1] "order" "n"    

# It is on the structural branch, so there is no block.
try(term_matrix(b))
#> Error : Can't find method for `term_matrix(<modelterms7::RegimeTerm>)`.
```
