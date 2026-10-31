# Shape Constraints for a Smooth

Construct the shapes that
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
imposes on a smooth. `shape_increasing()` and `shape_decreasing()` give
a monotone curve, with a slope of at least `min_slope` in absolute
value; `shape_convex()` and `shape_concave()` give a curve of constant
sign of curvature.

## Usage

``` r
shape_increasing(link = linkfunctions7::log_link(), min_slope = 0)

shape_decreasing(link = linkfunctions7::log_link(), min_slope = 0)

shape_convex(link = linkfunctions7::log_link())

shape_concave(link = linkfunctions7::log_link())
```

## Arguments

- link:

  A linkfunctions7 link onto \\(0, \infty)\\, by default
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).
  [`linkfunctions7::softplus_link()`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.html)
  is the other shipped choice.

- min_slope:

  A single nonnegative number, the smallest slope in absolute value.

## Value

An S7 object of class `ShapeMonotone` (`shape_increasing()`,
`shape_decreasing()`) or `ShapeCurvature` (`shape_convex()`,
`shape_concave()`), both inheriting from
[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md).

## Details

A monotone shape parametrizes the slopes of the B-spline coefficients on
the Greville abscissae, \\s_j = \pm(\delta + g^{-1}(\gamma_j))\\ with
\\\delta\\ = `min_slope`. A curvature shape keeps the slope at the left
end free and parametrizes the second divided differences, \\\pm
g^{-1}(\gamma_j)\\. See
[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md)
for the map.

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

# A link onto the whole real line is rejected.
try(shape_concave(link = linkfunctions7::identity_link()))
#> Error : the link 'identity' maps onto (-Inf, Inf); a shape needs a link onto (0, Inf),
#>   such as log_link() or softplus_link().
```
