# Shape-Constrained Smooth of One Covariate

Constructs a smooth function of one covariate constrained to a shape:
monotone
([`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md),
[`shape_decreasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)),
of constant curvature
([`shape_convex()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md),
[`shape_concave()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)),
both
([`shape_increasing_convex()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
and the other three combinations), or unimodal or sigmoid
([`shape_unimodal()`](https://statmodels7.github.io/modelterms7/reference/shape_located.md),
[`shape_sigmoid()`](https://statmodels7.github.io/modelterms7/reference/shape_located.md)).
The term enters the formula of any distribution parameter, like
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md).

## Usage

``` r
scs(
  x,
  shape = shape_increasing(),
  smoother = basis7::bspline_smooth(),
  by = NULL,
  by_hyper = c("shared", "level"),
  hyper = NULL,
  id = NULL,
  label = NULL,
  sparse = NULL
)
```

## Arguments

- x:

  The covariate, a bare name or an expression in the data.

- shape:

  A
  [shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md),
  by default
  [`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md).

- smoother:

  A basis7 smoother whose basis is a B-spline:
  [`basis7::bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.html)
  (the default) or
  [`basis7::pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.html).
  Only its basis is used.

- by:

  An optional factor or numeric variable, given as a bare name or an
  expression in the data, as in
  [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md).

- by_hyper:

  With a factor `by`, whether the levels share one smoothing parameter
  (`"shared"`, the default) or each has its own (`"level"`). Under
  `"level"` a held value or an `id` is qualified by the level,
  `hyper = c(lambda.a = 2)`.

- hyper:

  `NULL` to estimate the smoothing parameter, or a value at which to
  hold it, as for
  [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md).

- id:

  `NULL`, or a string shared with other terms whose smoothing parameter
  is the same, as for
  [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md).

- label:

  The term's label; `NULL` builds `scs(<x>)` from the call.

- sparse:

  `NULL` or `FALSE`: the block is dense.

## Value

An unbuilt
[ScsTerm](https://statmodels7.github.io/modelterms7/reference/ScsTerm.md).

## The construction

The covariate is expanded in the B-spline basis of `smoother`, whose
coefficients are written as \\\beta = A\\b(\gamma)\\ through the map of
the shape (see
[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md)):
cumulative sums of slopes, or of second divided differences, on the
Greville abscissae of the basis, each constrained one through the
inverse of the shape's link. Every value of the free coefficients
\\\gamma\\ gives a curve of the shape, so the constraint holds exactly
at every iterate of a fit.

The level is removed against the intercept of the model: the first
column of \\A\\ is dropped and the remaining columns of \\X A\\ are
centered over the data on which the term is built, so the smooth
averages to zero there. The coefficients are therefore \\\gamma_2,
\dots, \gamma_K\\, named `<label>.g1`, `<label>.g2`, and so on.

A shape with a value at a point, `value = c(at = a, value = v)`, removes
the level by that value instead: the columns are taken relative to their
row at \\a\\, and the contribution is \\f(x) = v + (B(x) -
B(a))\\A\_{-1}\\b\_{-1}(\gamma)\\. A convex or concave shape with a
slope at a point, `slope = c(at = a, value = v)`, has one coefficient
fewer, because the slope at the left end is solved from \\f'(a) = v\\;
the line that this slope contributes is a fixed part of the
contribution. See
[`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
for when a value other than zero changes the fit.

## The penalty

The penalty is the quadratic penalty of first differences of the
constrained coordinates of \\\gamma\\, under one smoothing parameter
`lambda` selected with the rest of the model. For a convex or concave
smooth these are the second divided differences, and the slope at the
left end stays free; for a combination the anchored slope is constrained
and is the first coordinate of the differences. Under a strong penalty a
monotone smooth approaches a straight line, and a convex or concave one
the curve whose coefficients have constant second divided differences,
which is a parabola in the interior and has a smaller curvature in the
end panels of the basis, where the knots are repeated. The penalty of
the smoother is not used: the smoother supplies the basis only.

## The block

The contribution is nonlinear in \\\gamma\\, so the block is the
Jacobian \\X_c A\_{-1}\\\mathrm{diag}(b'(\gamma))\\ at the current
coefficients, refreshed by
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
and
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
gives the contribution itself. Because \\b\\ acts on one coordinate at a
time, the derivatives of the block are diagonal in \\\gamma\\.

## Prediction outside the data

Beyond the ends of the basis the smooth is continued by the straight
line with the slope at the end, which preserves the shape.

## Positive and bounded functions

The shapes do not include a positive or a bounded function. When a
parameter must be positive, or must lie between two values, the range is
set by the link of that parameter in the distribution, such as
[`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
or
[`linkfunctions7::bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.html),
and [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md)
gives the smooth on the scale of the link. The whole predictor then lies
in the range at every iterate. A positive term beside an intercept would
not be identified, because a constant moved between the intercept and
the term leaves the likelihood unchanged.

## Varying the shape by another variable

A factor `by` gives one constrained curve per level, each with the shape
of `shape`: the block of level \\l\\ is the block of one curve
multiplied by the indicator of that level, every level is centered over
all the rows, as for
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md), and
the level means come from the factor written in the formula. The
coefficients are named `<label>.<level>.g1`, `<label>.<level>.g2`, and
so on. A value or a slope at a point holds for every level. Under
`by_hyper = "level"` each level has its own smoothing parameter.

A numeric `by` \\z\\ gives a varying-coefficient term \\z\\f(x)\\, with
the shape on the coefficient \\f\\. The term is not centered: the level
of \\f\\ is the coefficient of \\z\\, a free and unpenalized coordinate
named `<label>.const`, so writing \\z\\ in the formula as well makes it
aliased. With a value at a point the level is fixed by that value
instead, and `value = c(at = a, value = 0)` sets the effect of \\z\\ to
zero at \\a\\. The default label is `scs(x):z`.

Under a factor `by` the shape may also differ between the levels: its
declared settings (a minimum slope, a location, a peak or a valley, the
point and the value of a value or a slope) may be given by level, and
its estimated ones (a location that is not given, the left slope of a
curvature shape) may be one for every level or one per level, through
the arguments `mode_by`, `inflection_by` and `left_slope_by` of the
shape (see
[`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
and
[`shape_unimodal()`](https://statmodels7.github.io/modelterms7/reference/shape_located.md)).
A location chosen by level is searched by
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
one level at a time from a pilot fit of `s(x, by = g)`
([`term_variant_stages()`](https://statmodels7.github.io/modelterms7/reference/term_variant_stages.md)).

## Inference

Pointwise intervals from the curvature at the mode cover somewhat less
than their nominal level where the curve is flat (about 93 per cent for
95 per cent in simulation);
[`statmodels7::statmod_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/statmod_bootstrap.html)
gives a reading that does not rest on that curvature.

## References

Pya, N. and Wood, S. N. (2015). Shape constrained additive models.
*Statistics and Computing*, 25, 543–559.

## See also

[`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
and the other shapes,
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) for
the unconstrained smooth.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(100))
tm <- term_build(scs(x, shape_increasing(), basis7::bspline_smooth(k = 8)), dd)
term_coef_names(tm)
#> [1] "scs(x).g1" "scs(x).g2" "scs(x).g3" "scs(x).g4" "scs(x).g5" "scs(x).g6"
#> [7] "scs(x).g7"

# Every coefficient vector gives an increasing curve.
f <- term_value(tm, coef = rnorm(7))
all(diff(f[order(dd$x)]) >= -1e-12)
#> [1] TRUE
```
