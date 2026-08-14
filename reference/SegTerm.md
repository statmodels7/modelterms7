# S7 Class for Break-Point Terms

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a covariate whose effect changes at estimated break-points: a change
of slope
([`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md)), a
change of level
([`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)),
or both at the same points
([`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)).
The design block is the working one of the iteration that estimates the
break-points, and is recomputed by
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
as they move.

## Usage

``` r
SegTerm(
  label = character(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  kind = character(0),
  var = NULL,
  npsi = integer(0),
  linear = logical(0),
  subformulas = list(),
  spec = list()
)
```

## Arguments

- label:

  A character string prefixed to the coefficient names when non-empty.

- X:

  The design block, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- coef_names:

  The coefficient names, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- blueprint:

  The information needed to reproduce the mapping on new data, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- penalty:

  The penalty on the term's coefficients, or `NULL`.

- kind:

  Which of the three constructions.

- var:

  The covariate expression.

- npsi:

  The number of break-points.

- linear:

  Whether the block carries the linear effect.

- subformulas:

  A named list of one-sided formulas, one per developed parameter.

- spec:

  The resolved construction settings.

## Value

An object of class `SegTerm`.

## See also

[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)

## Examples

``` r
S7::S7_inherits(seg(x), SegTerm)
#> [1] TRUE
```
