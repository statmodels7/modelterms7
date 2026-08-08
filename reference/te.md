# Penalized Smooth of Several Covariates

A tensor-product smooth: a basis7 basis in each covariate, combined by
[`tensor_basis`](https://statmodels7.github.io/basis7/reference/tensor_basis.html),
with a roughness penalty on the product coefficients.

## Usage

``` r
te(..., by = NULL, k = 5, degree = 3, bases = NULL, label = NULL)
```

## Arguments

- ...:

  The covariates, expressions evaluated in the data, at least two of
  them.

- by:

  An optional factor or numeric variable, as in
  [`s`](https://statmodels7.github.io/modelterms7/reference/s.md).

- k:

  The basis dimension per margin, recycled to the number of covariates.
  Defaults to 5.

- degree:

  The spline degree per margin, recycled. Defaults to 3.

- bases:

  An optional list of basis7 bases, one per covariate, used in place of
  the default B-splines.

- label:

  A single non-empty string prefixed to the coefficient names. Defaults
  to a name built from the covariates.

## Value

An object of class
[`SmoothTerm`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The block is the tensor basis evaluated at the covariates, and the
penalty is the sum of the marginal roughness penalties carried into the
product, \\P = \sum_v I \otimes \cdots \otimes P_v \otimes \cdots
\otimes I\\, which penalizes curvature in each direction. One smoothing
parameter governs the sum, so the smoothing is isotropic across the
margins after each has been scaled to a common size; a separate
parameter per margin needs a penalty that is a sum of quadratics with
its own coefficient on each, which penalties7 does not provide.

The marginal bases are not reparametrized, so unlike
[`s`](https://statmodels7.github.io/modelterms7/reference/s.md) the
linear effects are not separated out: the null space of the tensor
penalty contains the constant and the marginal linear terms, and a model
carrying an intercept should constrain the smooth or accept that the
constant is shared.

## Examples

``` r
set.seed(2)
dd <- data.frame(x = runif(120), z = runif(120))
dd$y <- dd$x * dd$z + rnorm(120, sd = 0.1)
built <- term_build(te(x, z, k = 4), dd)
term_npar(built)
#> [1] 16
```
