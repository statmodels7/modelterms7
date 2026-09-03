# S7 Class for Structural Terms

The branch of
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
whose terms rewrite the likelihood instead of adding design columns.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
drives one distribution parameter by a score-driven recursion, so its
contribution to the predictor is a state rather than a product
\\X\beta\\;
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
and the marginal break-point terms replace the log-likelihood outright
with one mixed over latent states. A term on this branch has no block,
no coefficients and no
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
method.

## Usage

``` r
structural_term(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  ids = character(0),
  blueprint = list()
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

  A named list of everything needed to reproduce the term on new rows:
  the design's state, the levels a factor had, whatever the subclass's
  own filter or recursion needs. Its contents are that subclass's
  business and nothing outside the term reads them. Empty in a
  specification and filled by
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
  which is what
  [`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
  reads on this branch, a structural term having no coefficient names to
  record being built in.

## Value

Nothing: the class is abstract and cannot be instantiated.
`structural_term()` throws. As a type it is the parent of
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
and the marginal break-point terms, and carries
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)'s
six properties together with `blueprint`, the branch's record of having
been built, empty in a specification.

## Parameters, not coefficients

What such a term estimates are its **own parameters**, named by
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
and each riding a chart named by
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md),
so a persistence stays inside \\(-1, 1)\\ and a loading stays positive
whatever an optimizer proposes.
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
counts those parameters,
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
says where they start, and a
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
entry indexes into them, where an additive term's entry indexes into
columns.

The class adds no properties of its own: everything a built structural
term records goes in the root's `label` and in the subclass's own slots.

## The two shapes

A **filter** reports a predictor and its exact derivative:
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
returns the level at each observation together with the Jacobian in the
term's parameters, propagated beside the state because the recursion is
the only place it can be computed.
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md),
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
and
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
carry the first three orders through the same recursion.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
has this shape.

A **likelihood** term reports no predictor at all.
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
returns its own contribution,
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
the smoothed latent states and
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
the observed information of the mixture.
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
and the marginal break-point terms have this shape.

A fitting layer tells them apart by which of the two generics answers,
and at most one structural term is allowed per model formula.

## See also

[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for the other branch,
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
and
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
for the two shapes,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
and
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for what such a term estimates.

## Examples

``` r
# gas() and regime() are on this branch; every additive term is not.
vapply(list(gas(p = 1, q = 1), regime(k = 2), s(x, k = 5), ridge(~ x)),
       function(t) S7::S7_inherits(t, structural_term), logical(1))
#> [1]  TRUE  TRUE FALSE FALSE

# A structural term names its own parameters instead of coefficients.
g <- gas(p = 1, q = 2)
term_params(g)
#> [1] "omega"  "alpha1" "pacf1"  "pacf2" 
vapply(term_links(g), function(l) l@link_name, character(1))
#>      omega     alpha1      pacf1      pacf2 
#> "identity"      "log"   "rhobit"   "rhobit" 

# And has no design block at all.
try(term_matrix(g))
#> Error : Can't find method for `term_matrix(<modelterms7::GasTerm>)`.

# Abstract: there is nothing to construct.
try(structural_term())
#> Error in S7::new_object(S7::S7_object(), label = label, hyper = hyper,  : 
#>   Can't construct an object from abstract class <structural_term>
```
