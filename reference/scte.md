# Shape-Constrained Tensor Product Smooth of Two Covariates

Constructs a smooth function of two covariates that is monotone in one
of them for every value of the other, or monotone in both. The term
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
  label = NULL
)
```

## Arguments

- x, z:

  The covariates, bare names or expressions in the data.

- shapes:

  A list of two elements, for `x` and for `z`: each is
  [`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md),
  [`shape_decreasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
  or `NULL` for a free margin, and at least one is a shape. A shape here
  takes neither `min_slope` nor a value or a slope at a point.

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

- label:

  The term's label; `NULL` builds `scte(<x>, <z>)` from the call.

## Value

An unbuilt
[ScteTerm](https://statmodels7.github.io/modelterms7/reference/ScteTerm.md).

## The construction

Each covariate is expanded in the B-spline basis of its smoother, and
the surface is \\f(x, z) = B_x(x)^\top \Theta B_z(z)\\ with \\\Theta\\
the \\K_x \times K_z\\ matrix of the tensor-product coefficients. A
margin with a shape writes its direction of \\\Theta\\ as cumulative
sums of increments on the Greville abscissae of its basis, as
[`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
does for one covariate, each increment constrained through the inverse
of the log link. With the shape on \\x\\ alone, \\\Theta = A_x B\\:
every column of \\\Theta\\, a function of \\x\\ at one coefficient of
\\z\\, is monotone, and the first row of \\B\\ is a free function of
\\z\\ whose coefficients sum to zero, the constant being left to the
intercept of the model. With shapes on both margins, \\\Theta = A_x B
A_z^\top\\ and every element of \\B\\ except the level \\b\_{11}\\ is
constrained. A decreasing margin is the increasing one in the reflected
covariate, which reverses the order of its basis.

## The interaction under two shapes

With shapes on both margins every increment is constrained, the cross
increments \\b\_{jk}\\, \\j, k \ge 2\\, included. The surface is then
monotone in both covariates, and its interaction has a fixed sign: the
mixed differences of \\\Theta\\ are nonnegative when both margins are
increasing or both decreasing, and nonpositive when one is increasing
and the other decreasing. A surface that is increasing in both
covariates with a negative interaction, such as \\x + z - xz/2\\, is not
of this form, and the fit is biased for it; with the shape on one margin
the interaction is free. The construction is that of Pya and Wood
(2015).

## The penalty

Two smoothing parameters, `lambda1` and `lambda2`, are selected with the
rest of the model. `lambda1` penalizes the first differences of the
constrained coordinates of \\\gamma\\ along \\x\\, and `lambda2` those
along \\z\\; the free row or column has the second differences of its
coefficients along its own direction under the smoothing parameter of
that direction.

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
```
