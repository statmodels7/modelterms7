# S7 Class for a Smooth Random Effect

The subclass of
[`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md)
whose within-group design is a smooth: one random curve or surface per
level of the grouping variable.
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
constructs it when the left side of the bar is an
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) or a
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)
call, as in `random(~ s(x) | g)`.

## Usage

``` r
RandomSmoothTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  ids = character(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  formula = NULL,
  correlated = logical(0),
  distrib = NULL,
  tag = character(0),
  smooth = NULL,
  intercept = logical(0)
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

  The hyperparameters of the term's penalty that the caller **held**, as
  a named list keyed by the penalty's names. Empty, the default, means
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

- tag:

  The covariance label from the middle of two bars, or `NA_character_`.

- smooth:

  The unbuilt
  [`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
  of the left side of the bar.

- intercept:

  A single logical: whether each group has a random level.

## Value

An S7 object of class `RandomSmoothTerm`, inheriting from
[`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md),
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md).

## The two properties of its own

`smooth` is the unbuilt
[`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
read from the left side of the bar, which carries the basis, the
roughness penalty and the null space through its basis7 smoother.
`intercept` records whether each group also has a random level: `TRUE`
for `~ s(x) | g` and `FALSE` for `~ 0 + s(x) | g`.

Every other property, and every generic that
[`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md)
implements, has the meaning documented there: the block is sparse and
ordered group by group, a level that the build did not see is predicted
at a zero effect on request, and
[`term_within()`](https://statmodels7.github.io/modelterms7/reference/term_within.md)
returns the design that one group's effect multiplies.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
the constructor, whose section on smooth random effects gives the model;
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) for
the smooth itself.

## Examples

``` r
tm <- random(~ s(x, smoother = basis7::bspline_smooth(k = 6)) | g)
S7::S7_inherits(tm, RandomSmoothTerm)
#> [1] TRUE
S7::S7_inherits(tm, RandomTerm)
#> [1] TRUE
tm@intercept
#> [1] TRUE
```
