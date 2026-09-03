# S7 Class for Score-Driven Dynamics

The subclass of
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
holding a generalized autoregressive score component: a time-varying
level driven by the score of the observation density, added to the
predictor of one distribution parameter.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
constructs it. Its contribution is a state, so it implements
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
and has no
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
method at all.

## Usage

``` r
GasTerm(
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
  submodels = list()
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

  The number of score lags, an integer of at least 0.

- q:

  The number of autoregressive lags, an integer of at least 0.

- by:

  An optional grouping expression; each group is filtered independently.
  `NULL` for one series.

- time:

  An optional ordering expression. `NULL` for row order.

- links:

  A named list of linkfunctions7 links overriding the defaults, empty
  where none was given.

- submodels:

  A named list of one-sided formulas, one per parameter developed over
  covariates. Empty where none is.

## Value

An S7 object of class `GasTerm`, inheriting from
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the seven properties above beside
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)'s
six.

## The seven properties of its own

`p` is the number of score lags and `q` the number of autoregressive
ones, which together fix the parameter count at \\1 + p + q\\ before any
subformula.

`by` and `time` are the grouping and ordering expressions as written,
kept unevaluated; `NULL` means one series in row order.

`links` holds whatever the caller overrode, empty where the defaults
stand: the identity on the level, the log on each loading, the rhobit on
each partial autocorrelation. `submodels` holds one right-hand side per
parameter developed over covariates.

`blueprint` is filled by
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
and carries the row order within each group, the built sub-terms of each
subformula and their designs. The class overrides the branch's
`blueprint` property because a structural term has no design block to
hang one on.

## What the class is for

The recursion, its exact Jacobian, the reverse pass, the curvature and
the contracted third derivative are all methods on it, so a fitting
layer can estimate the term's parameters beside the coefficients of
every equation and read the joint observed information.

## See also

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
the constructor;
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
for the recursion;
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
for the quantities a fitted one reports;
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
for the other dynamic term.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:60, y = c(rnorm(30), rnorm(30, 3)))

tm <- gas(p = 1, q = 2, time = t)
S7::S7_inherits(tm, GasTerm)
#> [1] TRUE
c(p = tm@p, q = tm@q)
#> p q 
#> 1 2 

# 1 + p + q parameters, and one chart each.
term_params(tm)
#> [1] "omega"  "alpha1" "pacf1"  "pacf2" 
vapply(term_links(tm), function(l) l@link_name, character(1))
#>      omega     alpha1      pacf1      pacf2 
#> "identity"      "log"   "rhobit"   "rhobit" 

# The build resolves the ordering; there is no block to read.
b <- term_build(tm, dd)
names(b@blueprint)
#> [1] "order"  "n"      "levels" "times"  "group"  "by"     "time"  
try(term_matrix(b))
#> Error : Can't find method for `term_matrix(<modelterms7::GasTerm>)`.
```
