# Unimodal and Sigmoid Shapes

Construct the shapes of
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
that have a location: `shape_unimodal()` rises to a mode and falls after
it (or, with `valley = TRUE`, falls to a minimum and rises after it);
`shape_sigmoid()` increases (or, with `decreasing = TRUE`, decreases)
with a slope that grows up to an inflection and shrinks after it.

## Usage

``` r
shape_unimodal(
  link = linkfunctions7::log_link(),
  mode = NULL,
  index = NULL,
  valley = FALSE,
  value = NULL,
  mode_by = c("level", "shared")
)

shape_sigmoid(
  link = linkfunctions7::log_link(),
  inflection = NULL,
  index = NULL,
  decreasing = FALSE,
  value = NULL,
  inflection_by = c("level", "shared")
)
```

## Arguments

- link:

  A linkfunctions7 link onto \\(0, \infty)\\.

- mode, inflection:

  `NULL`, a single number (the value of the covariate at which the mode
  or the inflection is placed), or a vector of such numbers named by the
  levels of a factor `by`.

- index:

  `NULL`, the index \\M\\ itself, or a vector of indices named by level.
  A fitting layer sets it on the versions it tries. It takes precedence
  over `mode` and `inflection`.

- valley:

  `FALSE` for a maximum, `TRUE` for a minimum, or a logical vector named
  by level.

- value:

  `NULL`, two finite numbers `c(at = a, value = v)`, or a list
  `list(at = a, value = v)` whose elements are one number or a vector
  named by level: the smooth is \\v\\ at \\a\\, in place of the
  centering (see
  [`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
  for when \\v\\ changes the fit).

- mode_by, inflection_by:

  Under a factor `by`, whether a location that is not given is chosen
  for each level (`"level"`, the default) or once for every level
  (`"shared"`). `"shared"` without a factor `by` is an error.

- decreasing:

  `FALSE` for an increasing sigmoid, `TRUE` for a decreasing one, or a
  logical vector named by level.

## Value

An S7 object of class `ShapeUnimodal` or `ShapeSigmoid`, both inheriting
from `ShapeLocated` and
[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md).

## Where the location sits

The location is an index \\M\\ of the Greville abscissae \\\xi_1 \<
\dots \< \xi_K\\ of the basis. For `shape_unimodal()` the slopes of the
coefficients are \\\pm g^{-1}(\gamma_j)\\, of one sign up to \\M\\ and
of the other after it, so the coefficients, and with them the curve,
have a single extremum; \\M = 1\\ and \\M = K\\ are monotone curves. For
`shape_sigmoid()` every slope has the sign of the shape: the slope at
each end is constrained, and the second divided differences are summed
from each end toward the inflection, up to index \\M\\ from the left and
down to \\M + 1\\ from the right, so the slopes grow in absolute value
up to the inflection and shrink after it; \\M\\ runs from 2 to \\K\\.

The penalty of
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md) is
one random walk on \\\gamma\\ across the location, whose null space is
the same for every \\M\\, so fits at different \\M\\ are compared by
their criterion.

## How the location is chosen

With `mode` or `inflection` given, \\M\\ is the index of the Greville
abscissa nearest to it, and the model is fitted from two starts (the
least-squares lines of
[`shape_start()`](https://statmodels7.github.io/modelterms7/reference/shape_start.md)
and the slopes of a smoothing spline) whose better criterion is kept:
away from the mode of the data the first start is nearly flat and can
leave the fit on the flat curve. Left at `NULL`, the fitting layer
chooses \\M\\: a pilot fit with
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) in
place of the term gives the abscissa nearest to its extremum (unimodal)
or to its steepest point (sigmoid), \\M_0\\, and the model is fitted at
the indices \\M_0 - 2, \dots, M_0 + 2\\ that exist; the best criterion
decides (see
[`term_pilot()`](https://statmodels7.github.io/modelterms7/reference/term_pilot.md)
and
[`term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.md)).

## Under a factor `by`

With a factor `by` in
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md),
`mode`, `inflection`, `index`, `valley` and `decreasing` may be given by
level, as a vector named by the levels
(`valley = c(a = FALSE, b = TRUE)` gives a maximum in level `a` and a
minimum in level `b`), and the point and the value of `value` may be
given by level in the list form. A location that is not given is chosen
by level (`mode_by` or `inflection_by` at `"level"`, the default): the
fitting layer moves the index of one level at a time, from the pilot
indices of the levels. With `"shared"` one index holds for every level,
and the versions are taken around the median of the levels' pilot
indices.

## See also

[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md),
[`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md).

## Examples

``` r
shape_unimodal()
#> <shape: unimodal, link = log>
shape_unimodal(mode = 0.4)
#> <shape: unimodal, link = log, located at 0.4>
shape_unimodal(valley = TRUE)
#> <shape: unimodal (valley), link = log>
shape_sigmoid(inflection = 0.5)
#> <shape: sigmoid, link = log, located at 0.5>
shape_sigmoid(decreasing = TRUE)
#> <shape: sigmoid (decreasing), link = log>

# Under a factor by with levels a and b: a peak in a, a valley in b, and
# one location for both.
shape_unimodal(valley = c(a = FALSE, b = TRUE), mode_by = "shared")
#> <shape: unimodal (peak or valley by level), link = log, one location for every level>
```
