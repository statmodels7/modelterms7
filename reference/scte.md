# Shape-Constrained Tensor Product Smooth of Two Covariates

Constructs a smooth function of two covariates that has a shape in one
of them for every value of the other, or a shape in each: monotone,
convex or concave, or monotone and convex or concave together. The term
enters the formula of any distribution parameter, like
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md).

## Usage

``` r
scte(
  x,
  z,
  shapes = list(shape_increasing(), NULL),
  smooths = basis7::bspline_smooth(k = 6),
  by = NULL,
  by_hyper = c("shared", "level"),
  hyper = NULL,
  id = NULL,
  point = NULL,
  label = NULL
)
```

## Arguments

- x, z:

  The covariates, bare names or expressions in the data.

- shapes:

  A list of two elements, for `x` and for `z`: each is a monotone shape
  ([`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md),
  [`shape_decreasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)),
  a curvature shape
  ([`shape_convex()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md),
  [`shape_concave()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md),
  [`shape_increasing_convex()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
  and the other three combinations) or `NULL` for a free margin, and at
  least one is a shape. A monotone shape may set `min_slope`, a single
  value, and then the derivative in its covariate is at least
  `min_slope` in absolute value everywhere. A shape may set `value`, the
  value of the surface along the line at that value of its covariate
  (see Details); it takes neither a slope at a point nor values by
  level.

- smooths:

  A basis7 smoother whose basis is a B-spline, or a list of two, one per
  covariate:
  [`basis7::bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.html)
  (the default, with six basis functions) or
  [`basis7::pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.html).
  Only the bases are used.

- by:

  An optional factor or numeric variable, given as a bare name or an
  expression in the data, as in
  [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md).

- by_hyper:

  With a factor `by`, whether the levels share the two smoothing
  parameters (`"shared"`, the default) or each level has its own pair
  (`"level"`).

- hyper:

  `NULL` to estimate the smoothing parameters, or values at which to
  hold them, named `lambda1` and `lambda2` (with the level appended
  under `by_hyper = "level"`, as `lambda1.a`).

- id:

  `NULL`, or strings shared with other terms whose smoothing parameters
  are the same, named as `hyper`.

- point:

  `NULL`, or the value of the surface at one point, as
  `c(<x> = a, <z> = b, value = v)` with the covariates as written in the
  call, or unnamed in that order.

- label:

  The term's label; `NULL` builds `scte(<x>, <z>)` from the call.

## Value

An unbuilt
[ScteTerm](https://statmodels7.github.io/modelterms7/reference/ScteTerm.md).

## The construction

Each covariate is expanded in the B-spline basis of its smoother, and
the surface is \\f(x, z) = B_x(x)^\top \Theta B_z(z)\\ with \\\Theta\\
the \\K_x \times K_z\\ matrix of the tensor-product coefficients. A
margin with a shape writes its direction of \\\Theta\\ through the map
\\A\\ of that shape on the Greville abscissae of its basis, as
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
does for one covariate: the coordinates that the shape constrains (the
increments of a monotone shape, the second differences of a convex or
concave one) are nonnegative through the inverse of the log link, and
the others (the level, and the slope of a plain curvature) are free. A
minimum slope bounds the elements of \\B\\ on a coordinate that is
constrained in that direction and free in the other: the slope of a
column of \\\Theta\\ is a convex combination of them, and the elements
constrained in both directions add to it.

With the shape on \\x\\ alone, \\\Theta = A_x B\\. Every column of
\\\Theta\\ has the shape, and since the basis in \\z\\ is nonnegative
and sums to one, \\f(\cdot, z)\\ is a convex combination of the columns
and has the shape for every \\z\\. The rows of \\B\\ on the free
coordinates of the shape are free functions of \\z\\; the coefficients
of the level row sum to zero, the constant being left to the intercept
of the model.

With shapes on both margins, \\\Theta = A_x B A_z^\top\\. Every column
of \\\Theta\\ has the shape in \\x\\ when an element \\b\_{jk}\\ on a
coordinate \\j\\ constrained in \\x\\ has the sign of column \\k\\ of
\\A_z\\, and every row has the shape in \\z\\ when an element on a
coordinate constrained in \\z\\ has the sign of its column of \\A_x\\;
the free coordinates of a plain curvature are written as the two hat
functions of the ends of its interval, which are nonnegative. The shapes
hold for every coefficient vector. An element constrained in both
directions must then have both signs; the covariates are reflected
(which reverses a monotonicity and keeps a curvature) and the surface
changes sign where this makes the two signs agree, which it does for
every pair except a convex margin with a concave one. For that pair the
elements constrained in both directions are zero, and the curvature in
one covariate varies with the other only linearly. The elements free in
both directions sum to zero, the constant being left to the intercept.

## The interaction under two shapes

With shapes on both margins the cross elements \\b\_{jk}\\ are
constrained as well, so the interaction has a fixed sign. Under two
monotone shapes the mixed differences of \\\Theta\\ are nonnegative when
both margins are increasing or both decreasing, and nonpositive when one
is increasing and the other decreasing. A surface that is increasing in
both covariates with a negative interaction, such as \\x + z - xz/2\\,
is not of this form, and the fit is biased for it; with the shape on one
margin the interaction is free. The construction of two monotone margins
is that of Pya and Wood (2015).

## A value along a line or at a point

The shape of a margin may carry a value, `value = c(at = a, value = v)`
as in
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md),
and then \\f(a, z) = v\\ for every \\z\\ (or \\f(x, a) = v\\ for every
\\x\\ when the shape is on \\z\\). With \\r\\ the row of the basis of
\\x\\ at \\a\\ times the map of its shape, the condition is \\r B = c\\
for the coordinates \\c\\ of the constant \\v\\ in the other margin. In
each column of \\B\\ that has a free element, one free element (the
level where it enters) is written in terms of the others, and the
condition is exact. In a column without a free element \\c\\ is zero and
every element with \\r_j \ne 0\\ is zero, since all of them have the
sign of the shapes. When the other margin is free every column has a
free element, and the value removes the free function of \\z\\ at \\x =
a\\. When the other margin has a shape, the elements on the coordinates
of \\x\\ that the basis at \\a\\ touches are zero in the columns
constrained in \\z\\: the surface is flat in \\z\\ over that band, and
under two monotone shapes of the same direction over the whole side of
\\a\\ towards the start of the basis. At an end of the interval the band
is the first coefficient alone.

The argument `point` fixes the value of the surface at one point, \\f(a,
b) = v\\: the contribution is \\v + (W(x, z) - W(a, b))\\b\\, the number
of coefficients is unchanged, and under a numeric `by` the level of
\\f\\ is no longer a coefficient. With an intercept in the model a value
other than zero moves only the level that the intercept would carry. The
term takes one value: a line or a point. A point or a line outside the
interval of the basis is an error that names its bounds.

## The penalty

Two smoothing parameters, `lambda1` and `lambda2`, are selected with the
rest of the model. `lambda1` penalizes the first differences of the
constrained coordinates of \\\gamma\\ along \\x\\, and `lambda2` those
along \\z\\; the free row or column has the second differences of its
coefficients along its own direction under the smoothing parameter of
that direction. Under two monotone shapes the elements constrained in
\\x\\ only (the slopes in \\x\\ at the level of \\z\\), in \\z\\ only,
and in both (the interaction) are differenced only within their own
kind, so under a strong penalty the surface tends to \\a x + b z + c\\x
z\\ with \\c\\ of the sign of the interaction, as for
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md).
Under a curvature the differences run through the whole matrix.

## Varying the surface by another variable

A factor `by` gives one surface per level, each with the shapes of
`shapes`, and the level means come from the factor written in the
formula; the coefficients are named `<label>.<level>.g1`, and so on.
Under `by_hyper = "level"` each level has its own pair of smoothing
parameters, named `lambda1.<level>` and `lambda2.<level>`. A numeric
`by` \\w\\ gives the varying-coefficient term \\w\\f(x, z)\\: the level
of \\f\\ is then a free and unpenalized coefficient, `<label>.const`,
and the default label is `scte(x, z):w`.

## Prediction outside the data

Beyond the interval of a basis the surface is held at its value at the
boundary, which preserves the shapes.

## References

Pya, N. and Wood, S. N. (2015). Shape constrained additive models.
*Statistics and Computing*, 25, 543–559.

## See also

[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
for one covariate,
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) for
the unconstrained tensor product.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150), z = runif(150))
tm <- term_build(scte(x, z, shapes = list(shape_increasing(), NULL)), dd)
f <- term_value(tm, coef = rnorm(length(term_coef_names(tm))),
                newdata = expand.grid(x = seq(0, 1, 0.1), z = 0.3))
all(diff(f) >= -1e-12)
#> [1] TRUE

# Increasing in both covariates.
scte(x, z, shapes = list(shape_increasing(), shape_increasing()))
#> <ScteTerm> 'scte(x, z)' (specification; call term_build() with data); increasing in x, increasing in z

# Convex in x, increasing and concave in z.
scte(x, z, shapes = list(shape_convex(), shape_increasing_concave()))
#> <ScteTerm> 'scte(x, z)' (specification; call term_build() with data); convex in x, increasing concave in z

# Increasing in x with f(0, z) = 0 for every z, and a value at a point.
scte(x, z, shapes = list(shape_increasing(value = c(at = 0, value = 0)), NULL),
     smooths = basis7::bspline_smooth(k = 6, lower = 0, upper = 1))
#> <ScteTerm> 'scte(x, z)' (specification; call term_build() with data); increasing in x, free in z
scte(x, z, shapes = list(shape_increasing(), shape_convex()),
     point = c(x = 0.5, z = 0.5, value = 1))
#> <ScteTerm> 'scte(x, z)' (specification; call term_build() with data); increasing in x, convex in z
```
