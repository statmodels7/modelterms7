# Shape Constraints for a Smooth

Construct the shapes that
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
imposes on a smooth. `shape_increasing()` and `shape_decreasing()` give
a monotone curve, with a slope of at least `min_slope` in absolute
value; `shape_convex()` and `shape_concave()` give a curve of constant
sign of curvature; the four combinations `shape_increasing_convex()`,
`shape_increasing_concave()`, `shape_decreasing_convex()` and
`shape_decreasing_concave()` give a monotone curve of constant sign of
curvature.

## Usage

``` r
shape_increasing(
  link = linkfunctions7::log_link(),
  min_slope = 0,
  value = NULL
)

shape_decreasing(
  link = linkfunctions7::log_link(),
  min_slope = 0,
  value = NULL
)

shape_convex(
  link = linkfunctions7::log_link(),
  value = NULL,
  slope = NULL,
  left_slope_by = c("level", "shared")
)

shape_concave(
  link = linkfunctions7::log_link(),
  value = NULL,
  slope = NULL,
  left_slope_by = c("level", "shared")
)

shape_increasing_convex(
  link = linkfunctions7::log_link(),
  value = NULL,
  left_slope_by = c("level", "shared")
)

shape_increasing_concave(
  link = linkfunctions7::log_link(),
  value = NULL,
  left_slope_by = c("level", "shared")
)

shape_decreasing_convex(
  link = linkfunctions7::log_link(),
  value = NULL,
  left_slope_by = c("level", "shared")
)

shape_decreasing_concave(
  link = linkfunctions7::log_link(),
  value = NULL,
  left_slope_by = c("level", "shared")
)
```

## Arguments

- link:

  A linkfunctions7 link onto \\(0, \infty)\\, by default
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).
  [`linkfunctions7::softplus_link()`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.html)
  is the other shipped choice.

- min_slope:

  A single nonnegative number, the smallest slope in absolute value, or
  a vector of such numbers named by the levels of a factor `by`.

- value:

  `NULL`, two finite numbers `c(at = a, value = v)`, or a list
  `list(at = a, value = v)` whose elements are one number or a vector
  named by level: the smooth is \\v\\ at \\a\\. Unnamed, the first
  number is the point.

- slope:

  `NULL`, or the same forms as `value`: the derivative of the smooth is
  \\v\\ at \\a\\.

- left_slope_by:

  Under a factor `by`, whether the left slope of a curvature shape is a
  coefficient of each level (`"level"`, the default) or one coefficient
  for every level (`"shared"`). `"shared"` without a factor `by` is an
  error.

## Value

An S7 object of class `ShapeMonotone` (`shape_increasing()`,
`shape_decreasing()`) or `ShapeCurvature` (the other six), both
inheriting from
[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md).

## Details

A monotone shape parametrizes the slopes of the B-spline coefficients on
the Greville abscissae, \\s_j = \pm(\delta + g^{-1}(\gamma_j))\\ with
\\\delta\\ = `min_slope`. A curvature shape keeps the slope at the left
end free and parametrizes the second divided differences, \\\pm
g^{-1}(\gamma_j)\\. A combination constrains the slope as well, at the
end where the slope is smallest in absolute value: the left end for an
increasing convex or a decreasing concave curve, the right end for an
increasing concave or a decreasing convex one, so that every slope keeps
its sign. See
[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md)
for the map.

## A value or a slope at a point

With `value = c(at = a, value = v)` the smooth satisfies \\f(a) = v\\
exactly, and this constraint replaces the centering of
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md).
With `slope = c(at = a, value = v)`, which only `shape_convex()` and
`shape_concave()` accept, it satisfies \\f'(a) = v\\ exactly: the slope
at the left end, the one free slope of these two shapes, is then written
in terms of the second divided differences. The other shapes constrain
every slope, so a slope at a point cannot be imposed on them by this
elimination. The point \\a\\ must lie in the interval of the basis,
which the `lower` and `upper` arguments of the smoother can widen.

A value other than zero changes the fit only when no other term of the
equation absorbs a constant (an intercept, the effects of a factor, an
uncentered term), as in a formula that starts with `0 +`. Otherwise the
fit is the same for every \\v\\, and \\v\\ only moves the level between
the intercept and the smooth; with \\v = 0\\ the intercept is the linear
predictor at \\a\\.

## Under a factor `by`

With a factor `by` in
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md),
`min_slope` may be a vector named by the levels, and the point and the
value of `value` or `slope` may be given by level in the list form,
`list(at = 0, value = c(a = 0, b = 1))`. A curvature shape has one left
slope per level by default; with `left_slope_by = "shared"` that slope
(the free one of `shape_convex()` and `shape_concave()`, the anchored
one of a combination) is one coefficient for every level. The anchored
slope of a combination is penalized, so tying it is not available with
`by_hyper = "level"` in
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md),
which would penalize it once per level; a slope at a point removes the
left slope and cannot be combined with a tie.

## See also

[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md),
[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md).

## Examples

``` r
shape_increasing()
#> <shape: increasing, link = log>
shape_decreasing(min_slope = 0.1)
#> <shape: decreasing, link = log, min_slope = 0.1>
shape_convex(link = linkfunctions7::softplus_link())
#> <shape: convex, link = softplus(a=1)>
shape_increasing_concave()
#> <shape: increasing concave, link = log>
shape_increasing(value = c(at = 0, value = 0))
#> <shape: increasing, link = log, f(0) = 0>
shape_convex(slope = c(at = 0.5, value = 0))
#> <shape: convex, link = log, f'(0.5) = 0>

# Settings by level, for a factor by with levels a and b.
shape_increasing(min_slope = c(a = 0, b = 0.2))
#> <shape: increasing, link = log, min_slope = (a = 0.0, b = 0.2)>
shape_increasing(value = list(at = 0, value = c(a = 0, b = 1)))
#> <shape: increasing, link = log, f(0) = (a = 0, b = 1)>
shape_convex(left_slope_by = "shared")
#> <shape: convex, link = log, one left slope for every level>

# A link onto the whole real line is rejected.
try(shape_concave(link = linkfunctions7::identity_link()))
#> Error : the link 'identity' maps onto (-Inf, Inf); a shape needs a link onto (0, Inf),
#>   such as log_link() or softplus_link().
```
