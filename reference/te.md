# Penalized Smooth of Several Covariates

A tensor-product smooth: a basis7 basis in each covariate, combined by
[`tensor_basis`](https://statmodels7.github.io/basis7/reference/tensor_basis.html),
with a roughness penalty on the product coefficients.

## Usage

``` r
te(
  ...,
  by = NULL,
  k = 5,
  degree = 3,
  bases = NULL,
  anisotropic = TRUE,
  label = NULL
)
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

- anisotropic:

  Keep a smoothing parameter per margin? Defaults to `TRUE`.

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
penalty is built from the marginal roughness penalties carried into the
product, \\P_v = I \otimes \cdots \otimes P_v \otimes \cdots \otimes
I\\, each penalizing curvature in one direction.

With `anisotropic = TRUE`, the default, those components enter
[`additive_penalty`](https://statmodels7.github.io/penalties7/reference/additive_penalty.html)
and keep a smoothing parameter each, so the surface may be rough in one
direction and smooth in another – which is the usual reason for fitting
a tensor smooth rather than an isotropic one. With `anisotropic = FALSE`
they are summed first and one parameter governs the total, which costs
one hyperparameter instead of one per margin.

The marginal bases are not reparametrized, so unlike
[`s`](https://statmodels7.github.io/modelterms7/reference/s.md) the
linear effects are not separated out: the null space of the tensor
penalty contains the constant and the marginal linear terms, and a model
carrying an intercept should constrain the smooth or accept that the
constant is shared.

## References

Wood, S. N. (2006). Low-rank scale-invariant tensor product smooths for
generalized additive mixed models. *Biometrics* 62, 1025-1036.

Wood, S. N. (2017). *Generalized Additive Models: An Introduction with
R*, 2nd edition. Chapman and Hall/CRC.

## See also

[`s`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`random`](https://statmodels7.github.io/modelterms7/reference/random.md),
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md)

## Examples

``` r
set.seed(2)
dd <- data.frame(x = runif(120), z = runif(120))
dd$y <- dd$x * dd$z + rnorm(120, sd = 0.1)
built <- term_build(te(x, z, k = 4), dd)
term_npar(built)
#> [1] 16
```
