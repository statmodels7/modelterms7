# S7 Class for Marginal Break-Point Terms

The subclass of
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
holding break-points that vary by group as latent variables **integrated
out** of the likelihood.
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
construct it when given `marginal = TRUE` together with a
`psi ~ random(~ 1 | g)` subformula.

Its contribution is a likelihood, not a predictor, so it implements
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
and the prior over the positions is part of that likelihood itself.
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
declares nothing, and the prior's parameters are estimated by plain
maximum likelihood.

## Usage

``` r
MarginalBreakTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  ids = character(0),
  blueprint = list(),
  kind = character(0),
  var = NULL,
  npsi = integer(0),
  linear = logical(0),
  group = NULL,
  prior = NULL,
  spec = list()
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

- blueprint:

  A named list recording the number of rows, the covariate, the rows of
  each group and the group labels, empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- kind:

  One of `"jump"`, `"seg"` or `"jseg"`.

- var:

  The covariate expression, unevaluated.

- npsi:

  The number of break-points per group, an integer between 1 and 8 for
  `"jump"` and exactly 1 for the other two.

- linear:

  Whether the term carries the linear effect as its own parameter,
  `beta`. It is `FALSE` for `"jump"`.

- group:

  The grouping expression, from the break-point's
  [`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
  subformula.

- prior:

  The latent's distribution: `NULL` for the Gaussian, or a
  distributions7 object with its location held at zero.

- spec:

  A named list of the resolved construction settings.

## Value

An S7 object of class `MarginalBreakTerm`, inheriting from
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the properties listed above beside the inherited ones.

## Properties

The class defines the properties `kind`, `var`, `npsi`, `linear`,
`group`, `prior` and `spec`, and inherits the others from
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
among them `blueprint`.

`kind` is `"jump"`, `"seg"` or `"jseg"`; `var` is the covariate
expression; `npsi` is the number of break-points per group; `linear`
records whether the term carries the linear effect as a parameter of its
own. It is `TRUE` for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
unless `linear = FALSE` was given, and `FALSE` for
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md).

`group` is the grouping expression, taken from the break-point's
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
subformula. `prior` is the distribution of the latent positions: `NULL`
for the Gaussian, or a distributions7 object where `random(distrib = )`
names one. The location of that distribution must be fixed at zero,
since `psi1.mean` carries the position.

`spec` holds the resolved construction settings. `blueprint` is
inherited from
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
and holds the grouping and the covariate that
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
records.

## Parameters

The parameters are numbered, one set per break-point. A step term with
one break-point has `psi1.mean`, `psi1.sd` and `delta1`, where
`psi1.mean` is the mean of the positions in the population, `psi1.sd` is
their standard deviation between groups (estimated through a log link)
and `delta1` is the change of level. A
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
term has `beta` for the linear effect, `psi1.mean`, `psi1.sd` and
`gamma1` for the change of slope. A
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term has `beta`, `psi1.mean`, `psi1.sd`, `gamma1` and `delta1`.

## Number of break-points

A step term takes at most eight break-points per group. The forward
recursion of the step kind costs \\n K 2^K\\ operations for \\n\\
observations and \\K\\ break-points. The fitting layer reads a posterior
over the \\2^K\\ side patterns and evaluates the family once per
pattern, and the number of patterns therefore limits \\K\\. The
continuous kinds take one break-point, because a product quadrature over
several would need a node set whose size grows as a power of their
number.

## See also

[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
for the constructors;
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
for the likelihood;
[`term_latent()`](https://statmodels7.github.io/modelterms7/reference/term_latent.md)
for the posterior positions.

## Examples

``` r
set.seed(1)
dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)

tm <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
S7::S7_inherits(tm, MarginalBreakTerm)
#> [1] TRUE
c(kind = tm@kind, npsi = tm@npsi, linear = tm@linear)
#>    kind    npsi  linear 
#>  "jump"     "1" "FALSE" 

# Numbered parameters: the prior's location and scale, and the change.
term_params(tm)
#> [1] "psi1.mean" "psi1.sd"   "delta1"   
vapply(term_links(tm), function(l) l@link_name, character(1))
#>  psi1.mean    psi1.sd     delta1 
#> "identity"      "log" "identity" 

# The prior is part of the likelihood, so nothing is declared penalized.
length(term_penalties(tm))
#> [1] 0

# A continuous kind adds the linear effect and the change of slope.
term_params(term_build(seg(x, psi ~ random(~ 1 | id), marginal = TRUE), dd))
#> [1] "beta"      "psi1.mean" "psi1.sd"   "gamma1"   
```
