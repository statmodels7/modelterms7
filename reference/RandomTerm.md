# S7 Class for Grouped Random-Effect Terms

The subclass of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
holding grouped coefficients with a distribution attached to them: the
within-group design interacted with the group indicators, one
coefficient per group and per within-group column, and a penalty
carrying that distribution's negative log-density.
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
constructs it.

## Usage

``` r
RandomTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  formula = NULL,
  correlated = logical(0),
  distrib = NULL
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

- X:

  The design block, one row per observation: a numeric matrix or a
  two-dimensional Matrix. Empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- coef_names:

  The block's coefficient names, one per column, prefixed by `label`
  when there is one. Filled by
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- blueprint:

  A named list of everything needed to reproduce the mapping on new
  rows. Its contents are the subclass's business; nothing outside the
  term reads them. Empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- penalty:

  A penalties7 penalty on the block's coefficients, or `NULL` when the
  term is unpenalized.

- formula:

  The bar formula, `~ 1 | g` or `~ x | g`, with the within-group design
  on the left and the grouping variable on the right.

- correlated:

  A single logical: whether the default Gaussian lets the within-group
  effects correlate. Read only when `distrib` is `NULL`.

- distrib:

  The effects' distribution, a distributions7 object or a list of them
  with one per within-group column, or `NULL` for the default Gaussian.

## Value

An S7 object of class `RandomTerm`, inheriting from
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the three properties above beside the ten they supply.

## The three properties of its own

`formula` is the bar formula as given, `~ 1 | g` or `~ x | g`, kept with
its environment. `correlated` says whether the **default** Gaussian lets
the within-group effects depend on each other; it is read only where
`distrib` is `NULL`, the two saying the same thing.

`distrib` is the effects' distribution as supplied, or `NULL` for the
default. What the build turns it into is a penalties7 penalty, read
through
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
or
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
so the hyperparameter names and their bounds all come from the
distribution.

## The block is sparse by construction

A row belongs to one group, so the block has a density of \\1/m\\ and is
always a `dgCMatrix`.
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
takes no `sparse` argument for that reason, and passing one is an error.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
the constructor;
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the entries a built one declares;
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what a fitted random effect spends.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(9), g = factor(rep(c("a", "b", "c"), 3)))

tm <- random(~ x | g)
S7::S7_inherits(tm, RandomTerm)
#> [1] TRUE
tm@formula
#> ~x | g
#> <environment: 0x55a7fcf95830>
tm@correlated
#> [1] TRUE

# The block is one diagonal block per group and is always sparse.
b <- term_build(tm, dd)
class(term_matrix(b))
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
term_coef_names(b)
#> [1] "random.a.(Intercept)" "random.a.x"           "random.b.(Intercept)"
#> [4] "random.b.x"           "random.c.(Intercept)" "random.c.x"          

# The hyperparameters are the effects' distribution's own.
term_penalty(b)@params
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
```
