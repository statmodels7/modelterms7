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

  A named list of the resolved grouping and interval structure, empty
  until
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
  `TRUE` for `"seg"` and `"jseg"`.

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
with the eight properties above beside
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)'s
six.

## The eight properties of its own

`kind` is `"jump"`, `"seg"` or `"jseg"`; `var` the covariate expression;
`npsi` the number of break-points per group; `linear` whether the term
carries the linear effect as a parameter of its own, which `seg` and
`jseg` do.

`group` is the grouping expression, taken from the break-point's
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
subformula. `prior` is the latent's distribution: `NULL` for the
Gaussian, or a distributions7 object where `random(distrib = )` named
one, and its location must be fixed at zero, `m1` carrying the position.

`spec` holds the resolved construction settings and `blueprint` the
grouping and the interval structure
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
worked out.

## The parameters

They are numbered, one set per break-point: `m1`, `tau1`, `delta1` for a
one-break-point step term, with `m` the prior's location, `tau` its
scale on a log chart, and `delta` the change of level. A continuous kind
adds `beta` for the linear effect and `gamma1` for the change of slope.

## What it costs

At most **eight** break-points. The forward recursion of the step kind
costs \\n K 2^K\\ and stays cheap well past that; what does not is the
fitting layer, which reads a posterior over the \\2^K\\ side patterns
and evaluates the family once per pattern. The continuous kinds stay at
one break-point, a product quadrature over more being far dearer still.

## See also

[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
for the constructors;
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
for the likelihood;
[`term_latent()`](https://statmodels7.github.io/modelterms7/reference/term_latent.md)
for the posterior positions a reader wants.

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
#> [1] "m1"     "tau1"   "delta1"
vapply(term_links(tm), function(l) l@link_name, character(1))
#>         m1       tau1     delta1 
#> "identity"      "log" "identity" 

# The prior is part of the likelihood, so nothing is declared penalized.
length(term_penalties(tm))
#> [1] 0

# A continuous kind adds the linear effect and the change of slope.
term_params(term_build(seg(x, psi ~ random(~ 1 | id), marginal = TRUE), dd))
#> [1] "beta"   "m1"     "tau1"   "gamma1"
```
