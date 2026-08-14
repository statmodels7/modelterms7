# S7 Base Class for Model Terms

The abstract root of the term hierarchy. A model term records what a
formula names: the recipe that turns a data frame into a contribution to
the model, together with the metadata a fit reads. A term as written in
a formula is a specification;
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
turns it into a built term carrying its design block.

The formula interpreter
([`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md))
recognizes a term by this class: any call in a formula that evaluates to
an object inheriting from `model_term` is treated as a term, so a term
class defined outside the package works in a formula without
registration.

## Usage

``` r
model_term(
  label = character(0),
  hyper = list(),
  grid = list(),
  min_ratio = numeric(0)
)
```

## Arguments

- label:

  A character string prefixed to the term's coefficient names when
  non-empty.

- hyper:

  The hyperparameters of the term's penalty that the caller HELD, as a
  named list. Empty, the default, means every one of them is estimated.
  See
  [`term_hyper`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md).

- grid:

  How many values a path visits for each of the term's hyperparameters,
  as a named list. Empty, the default, leaves it to the criterion. See
  [`term_grid`](https://statmodels7.github.io/modelterms7/reference/term_grid.md).

- min_ratio:

  How far down the path over the size of the kink reaches, as a fraction
  of the value that empties the block, or `numeric(0)` for the
  criterion's own. See
  [`term_path_min`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md).

## Value

An object inheriting from class `model_term`.

## See also

[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md),
[`structural_term`](https://statmodels7.github.io/modelterms7/reference/structural_term.md),
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)

## Examples

``` r
S7::S7_inherits(linpar(~1), model_term)
#> [1] TRUE
```
