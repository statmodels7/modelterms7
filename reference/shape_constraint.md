# S7 Class for Shape Constraints

`shape_constraint` is the parent class of the shapes that
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
imposes on a smooth:
[`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md),
[`shape_decreasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md),
[`shape_convex()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
and
[`shape_concave()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md).
A shape writes the B-spline coefficients of the smooth as \\\beta =
A\\b(\gamma)\\, where \\A\\ is a fixed matrix of cumulative sums on the
Greville abscissae of the basis and \\b\\ applies the inverse of a link
to the constrained coordinates one at a time, so every free vector
\\\gamma\\ gives a curve of that shape.

## Usage

``` r
shape_constraint(shape_name = character(0), link = NULL, shape_params = list())

ShapeMonotone(shape_name = character(0), link = NULL, shape_params = list())

ShapeCurvature(shape_name = character(0), link = NULL, shape_params = list())
```

## Arguments

- shape_name:

  A single string naming the shape.

- link:

  A linkfunctions7 link onto \\(0, \infty)\\.

- shape_params:

  A named list of the shape's own settings.

## Value

An S7 object of class `shape_constraint`.

## The map

With \\\xi_1 \< \dots \< \xi_K\\ the Greville abscissae of a B-spline
basis and \\h_j = \xi\_{j+1} - \xi_j\\, a monotone shape parametrizes
the slopes \\s_j = (\beta_j - \beta\_{j-1}) / h\_{j-1}\\ and a convex or
concave shape the second divided differences of the coefficients. A
B-spline whose coefficients have nonnegative slopes is nondecreasing,
and one whose slopes on the Greville abscissae are nondecreasing is
convex, which holds for any knot placement, including the knots repeated
at the ends that basis7 uses.

On a constrained coordinate \\b_j = \ell_j + g^{-1}(\gamma_j)\\, with
\\g\\ the link of the shape (a link onto \\(0, \infty)\\) and \\\ell_j\\
a lower bound (zero, or the minimum slope). On a free coordinate \\b_j =
\gamma_j\\. The derivatives of the map are therefore diagonal:
\\\partial^k b_j / \partial \gamma_j^k = (g^{-1})^{(k)}(\gamma_j)\\ and
every mixed derivative is zero.

## A new shape

A user-defined shape is a subclass of `shape_constraint` with methods
for
[`shape_map()`](https://statmodels7.github.io/modelterms7/reference/shape_map.md)
and
[`shape_start()`](https://statmodels7.github.io/modelterms7/reference/shape_start.md).

## See also

[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md),
[`shape_map()`](https://statmodels7.github.io/modelterms7/reference/shape_map.md),
[`shape_start()`](https://statmodels7.github.io/modelterms7/reference/shape_start.md).

## Examples

``` r
sh <- shape_increasing()
S7::S7_inherits(sh, shape_constraint)
#> [1] TRUE
sh@shape_name
#> [1] "increasing"
```
