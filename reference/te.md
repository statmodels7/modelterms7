# Penalized Smooth of Several Covariates

A tensor-product smooth: a basis7 basis in each covariate, combined by
[`basis7::tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.html),
with a roughness penalty on the product coefficients. By default each
margin keeps a smoothing parameter of its own, so the surface may be
rough in one direction and smooth in another.

## Usage

``` r
te(
  ...,
  by = NULL,
  k = 5,
  degree = 3,
  bases = NULL,
  anisotropic = TRUE,
  label = NULL,
  lambda = NULL,
  sparse = NULL
)
```

## Arguments

- ...:

  The covariates, bare expressions evaluated in the data, at least two
  of them. One throws
  `"'te' needs at least two covariates; use s() for one."`.

- by:

  An optional factor or numeric variable, as in
  [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md),
  with the same two readings and the same sparsity rule.

- k:

  The basis dimension per margin, `5` by default, recycled to the number
  of covariates. As in
  [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) it
  must exceed `degree`.

- degree:

  The spline degree per margin, `3` by default, recycled.

- bases:

  An optional list of basis7 bases, one per covariate, used in place of
  the default B-splines.

- anisotropic:

  `TRUE`, the default, for one smoothing parameter per margin; `FALSE`
  for one over their sum. Anything that is not a single logical throws.

- label:

  A single non-empty string prefixed to the coefficient names. `NULL`,
  the default, builds one from the covariates: `te(x,z)`.

- lambda:

  The smoothing parameters, held at the values given and **estimated**
  when left `NULL`, which is the default. An anisotropic product carries
  one per margin, so a vector of that length, or a named one holding
  some of them.

- sparse:

  `TRUE`, `FALSE`, or `NULL` to settle it at build. Only a factor `by`
  admits `TRUE`. See
  [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md).

## Value

An unbuilt
[`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md):
a specification, with `X`, `coef_names`, `blueprint` and `penalty` empty
until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills them.

## The block and its penalty

The block is the tensor basis evaluated at the covariates, and the
penalty is built from the marginal roughness penalties carried into the
product, \\P_v = I \otimes \cdots \otimes P_v \otimes \cdots \otimes
I\\, each penalizing curvature in one direction.

With `anisotropic = TRUE`, the default, those components go to
[`penalties7::additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.html)
and keep one smoothing parameter each, named `lambda1`, `lambda2`, and
so on. That is the usual reason for fitting a tensor smooth. With
`anisotropic = FALSE` they are summed first and one `lambda` governs the
total, which costs one hyperparameter instead of one per margin.

The marginal bases are **not** reparametrized, so the marginal linear
effects are not separated out as
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md)
separates its one. They lie in the null space of the tensor penalty, and
a strongly penalized fit is shrunk toward no surface at all rather than
toward a plane.

## Centering

The tensor product of the marginal bases contains the constant, which
the penalty's null space contains as well, so beside an intercept the
block would be rank deficient by exactly one and the penalty would not
cover the deficiency.

The block therefore carries the sum-to-zero constraint over the observed
covariates
([`basis7::constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.html)).
The term has **one column fewer** than the product of its marginal
dimensions, so `te(x, z, k = 4)` gives 15 and not 16; every column sums
to zero over the data it was built on, to machine precision; and the
penalty follows by congruence with its rank unchanged, the direction
removed having been one of its null directions.

The transform is stored in the blueprint and reapplied by
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md),
as the Demmler-Reinsch transform of
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) is.

The level of the surface is then the model's intercept, so a formula
that removes it, `y ~ te(x, z) - 1`, fits a surface constrained to
average zero over the data.

## References

Wood, S. N. (2006). Low-rank scale-invariant tensor product smooths for
generalized additive mixed models. *Biometrics* 62, 1025–1036.

Wood, S. N. (2017). *Generalized Additive Models: An Introduction with
R*, 2nd edition. Chapman and Hall/CRC.

## See also

[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) for
one covariate,
[`basis7::tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.html)
for the product,
[`penalties7::additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.html)
for the anisotropic penalty.

## Examples

``` r
set.seed(2)
dd <- data.frame(x = runif(120), z = runif(120))
dd$y <- dd$x * dd$z + rnorm(120, sd = 0.1)

# Four by four margins give fifteen columns: the centering removes one.
b <- term_build(te(x, z, k = 4), dd)
c(npar = term_npar(b), product = 4 * 4)
#>    npar product 
#>      15      16 

# Every column sums to zero over the data, to machine precision.
max(abs(colSums(term_matrix(b))))
#> [1] 1.427244e-15

# Anisotropic by default: one smoothing parameter per margin.
term_penalty(b)@penalty_name
#> [1] "additive [2 components]"
term_penalty(b)@params
#> [1] "lambda1" "lambda2"
term_penalty(term_build(te(x, z, k = 4, anisotropic = FALSE), dd))@params
#> [1] "lambda"

# Holding both of them.
term_hyper(te(x, z, k = 4, lambda = c(1, 5)))
#> [[1]]
#> [[1]]$lambda1
#> [1] 1
#> 
#> [[1]]$lambda2
#> [1] 5
#> 
#> 

# The centering transform is reapplied, not recomputed.
max(abs(term_predict(b, dd[1:10, ]) - term_matrix(b)[1:10, ]))
#> [1] 4.163336e-17

# One covariate is s(), not te().
try(te(x, k = 4))
#> Error : 'te' needs at least two covariates; use s() for one.


# Fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(5)
  fd <- data.frame(a = runif(300, -2, 2), b = runif(300, -2, 2))
  fd$y <- sin(fd$a) * cos(fd$b) + rnorm(300, sd = 0.3)
  ft <- statmodels7::statmod(y ~ te(a, b, k = 5),
                             distributions7::gaussian1_distrib(), fd)
  # one smoothing parameter per margin, against a truth of sin(a) cos(b)
  round(sqrt(mean((fitted(ft) - sin(fd$a) * cos(fd$b))^2)), 3)
}
#> [1] 0.081
```
