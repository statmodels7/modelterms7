# The Starting Value of a Shape Constraint

Returns the starting value of \\b\\ for a shape: constant increments,
read from a straight line (monotone shapes) or a parabola (curvature
shapes) fitted by least squares to a working response.

## Usage

``` r
shape_start(shape, x, target, xi, ...)
```

## Arguments

- shape:

  A
  [shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md).

- x:

  The covariate.

- target:

  A working response on the scale of the linear predictor, of the length
  of `x`, or `NULL`, in which case the start is a flat curve at the
  floor.

- xi:

  The Greville abscissae of the basis.

- ...:

  Passed to methods.

## Value

A numeric vector of length `length(xi)`, the starting \\b\\ on the
natural scale (the level first).

## Details

A monotone shape starts from equal slopes, the slope of the
least-squares line when its sign agrees with the shape. A curvature
shape starts from the least-squares parabola: its slope at the left end
and its second derivative. A slope or a curvature of the wrong sign, or
below the lower bound, is replaced by a floor of one hundredth of the
scale of the response over the range of the covariate, so that \\g(b_j -
\ell_j)\\ is finite.

## See also

[`shape_map()`](https://statmodels7.github.io/modelterms7/reference/shape_map.md),
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md).

## Examples

``` r
set.seed(1)
x <- runif(50); y <- 2 * x + rnorm(50, sd = 0.1)
shape_start(shape_increasing(), x, y, xi = seq(0, 1, length.out = 6))
#> [1] 0.006271613 2.006778926 2.006778926 2.006778926 2.006778926 2.006778926
```
