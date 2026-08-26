# S7 Base Class for Model Terms

The abstract root of the term hierarchy. A term records what a formula
names: the recipe turning a data frame into a contribution to the model,
together with the metadata a fit reads. A term written in a formula is a
**specification**, carrying only what the call said;
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
turns it into a **built** term carrying the design block or the state
the recipe produces on given data.

[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
recognizes a term by this class. Any call on the right of a formula
whose value inherits from `model_term` becomes a term, so a term class
defined outside the package works in a formula the day it is written,
with nothing to register.

## Usage

``` r
model_term(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0)
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

## Value

Nothing: the class is abstract and cannot be instantiated. As a type it
is the parent of every term, with the six properties above.

## The two branches

[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
contributes a block of design columns \\X_j\beta_j\\ to a linear
predictor.
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
rewrites the likelihood instead:
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)'s
predictor is a recursion,
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)'s
contribution is a likelihood mixed over latent states. The branch
decides which generics a consumer may call, so
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
routes on it.

## What the root carries, and why

Beyond `label` the properties are all about a hyperparameter path.
**Which hyperparameters are estimated is said by the term and by nothing
else**: the term is where the penalty is named, so it is where a held
value, a grid size, a written-out set of values, the depth of the path
and the way several hyperparameters are combined all belong. An outer
criterion is put to every hyperparameter of a model, smooth ones
included, and carries none of this.

Each of `hyper`, `grid` and `values` is a named list keyed by the
penalty's own hyperparameter names, and each is empty by default.
`min_ratio` and `search` are single values rather than one per
hyperparameter, because only the path over the size of the kink uses a
ratio: a bounded hyperparameter is swept over its own interval and a
shape over a geometric grid above its lower bound.

## It cannot be constructed

`model_term` is abstract, as are `additive_term` and `structural_term`.
`model_term()` throws
`"Can't construct an object from abstract class <model_term>"`. Use it
as a parent when writing a term class, and as the test
`S7::S7_inherits(x, model_term)` when asking whether an object is a
term.

A class inheriting from it supplies
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
at the very least: the method registered here throws, naming the class
that did not implement it. Most other generics carry a usable default,
so a new term class starts from a working object and overrides what it
has reason to.

## See also

[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
for the two branches,
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
for turning a specification into a built term,
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
for how a formula is read into terms, and
[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
for validating one.

## Examples

``` r
# Every term inherits from this class, whichever branch it is on.
vapply(list(linpar(~ 1), ridge(~ x), s(x, k = 5), gas(p = 1, q = 1)),
       function(t) S7::S7_inherits(t, model_term), logical(1))
#> [1] TRUE TRUE TRUE TRUE

# The branch is what a consumer routes on.
c(additive = S7::S7_inherits(s(x, k = 5), additive_term),
  structural = S7::S7_inherits(s(x, k = 5), structural_term))
#>   additive structural 
#>       TRUE      FALSE 
c(additive = S7::S7_inherits(gas(p = 1, q = 1), additive_term),
  structural = S7::S7_inherits(gas(p = 1, q = 1), structural_term))
#>   additive structural 
#>      FALSE       TRUE 

# The path properties are empty until a caller sets one.
r <- ridge(~ x)
lengths(list(hyper = r@hyper, grid = r@grid, values = r@values,
             min_ratio = r@min_ratio, search = r@search))
#>     hyper      grid    values min_ratio    search 
#>         0         0         0         0         0 

# Holding lambda puts it in `hyper`, where term_hyper() reads it.
term_hyper(ridge(~ x, lambda = 0.5))
#> [[1]]
#> [[1]]$lambda
#> [1] 0.5
#> 
#> 

# Abstract: there is nothing to construct.
try(model_term())
#> Error in S7::new_object(S7::S7_object(), label = label, hyper = hyper,  : 
#>   Can't construct an object from abstract class <model_term>
```
