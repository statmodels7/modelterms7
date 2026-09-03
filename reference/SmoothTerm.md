# S7 Class for Smooth Terms

The subclass of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a penalized smooth: a basis7 expansion of one or more covariates
with a roughness penalty on its coefficients.
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md)
constructs it for one covariate and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) for
several, and the two differ in the basis they build and the penalty they
attach, not in the class.

## Usage

``` r
SmoothTerm(
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
  vars = list(),
  by = NULL,
  spec = list(),
  sparse = NULL
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

- vars:

  A list of the covariate expressions being smoothed.

- by:

  The expression the smooth varies with, or `NULL`.

- spec:

  A named list of construction settings: the basis, its dimension and
  degree, whether the linear part is carried separately, and for
  [`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)
  the `anisotropic` flag.

- sparse:

  `TRUE`, `FALSE` or `NULL` for the block's storage; only a factor `by`
  admits `TRUE`. See
  [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md).

## Value

An S7 object of class `SmoothTerm`, inheriting from
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the four properties above beside the ten they supply.

## The four properties of its own

`vars` holds the covariate expressions, one for
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
two or more for
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md),
kept unevaluated so that a build reads them in whatever data it is
given. `by` is the expression the smooth varies with, or `NULL`.

`spec` carries the construction settings the build reads: the basis or
bases, the dimension `k`, the degree, whether the linear part is
separated out, and for
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)
whether the penalty is anisotropic. What a build then computes from the
data goes into the blueprint instead: the Demmler-Reinsch transform, the
centering constraint, the `by` levels.

`sparse` is `NULL` until the build settles it. A smooth's block is
sparse only under a **factor** `by`, where each row sits in the block of
its own level: the basis itself is dense by construction, the
Demmler-Reinsch rotation making it so, and a numeric `by` merely
multiplies it.

## See also

[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md), the
two constructors;
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
for the roughness penalty;
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what a fitted smooth spends.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(80)), z = runif(80))

# Both constructors return this class.
c(s = S7::S7_inherits(s(x), SmoothTerm),
  te = S7::S7_inherits(te(x, z), SmoothTerm))
#>    s   te 
#> TRUE TRUE 

# The settings are on `spec`; what the data decide is in the blueprint.
tm <- s(x, k = 8)
names(tm@spec)
#> [1] "k"      "degree" "bases"  "linear"
names(term_build(tm, dd)@blueprint)
#> [1] "core"      "marg"      "spec"      "vars"      "by"        "by_levels"
#> [7] "sparse"    "nblock"   
```
