# The Component Matrices of a Built Smooth

Returns the roughness matrices of a built
[`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
as a list, one per smoothing parameter: one matrix for a single penalty,
several for an adaptive smoother or an anisotropic tensor product.

## Usage

``` r
.random_smooth_components(sb)
```

## Arguments

- sb:

  A built
  [`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
  without a `by` variable.

## Value

A list of symmetric numeric matrices of side
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

## Details

A smoother carrying a penalty factory has no roughness matrix, because
the factory declares its penalty over the penalized coordinates, and the
function signals an error for such a smoother.
