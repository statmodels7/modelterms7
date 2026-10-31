# The Map of a Shape Constraint

Returns the linear part of the map \\\beta = A\\b(\gamma)\\ of a shape
on a B-spline basis with Greville abscissae `xi`: the matrix \\A\\, the
positions of the constrained coordinates and their lower bounds.

## Usage

``` r
shape_map(shape, xi, ...)
```

## Arguments

- shape:

  A
  [shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md).

- xi:

  The Greville abscissae of the basis, a strictly increasing numeric
  vector of length \\K \ge 3\\.

- ...:

  Passed to methods.

## Value

A list with `A`, a \\K \times K\\ matrix; `cons`, the integer positions
of the constrained coordinates; and `lower`, a numeric vector of length
\\K\\, the bound \\\ell_j\\ added to \\g^{-1}(\gamma_j)\\ on a
constrained coordinate (zero elsewhere).

## Details

The first coordinate is always the level, \\\beta_1 = b_1\\, and the
first column of \\A\\ is a column of ones, because a B-spline basis sums
to one.
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
removes that column against the intercept of the model.

## See also

[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md),
[`shape_start()`](https://statmodels7.github.io/modelterms7/reference/shape_start.md).

## Examples

``` r
xi <- c(0, 0.2, 0.5, 0.8, 1)
m <- shape_map(shape_increasing(), xi)
m$cons
#> [1] 2 3 4 5

# Every positive b gives increasing coefficients.
b <- c(0.3, exp(rnorm(4)))
all(diff(drop(m$A %*% b)) > 0)
#> [1] TRUE
```
