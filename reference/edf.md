# Effective Degrees of Freedom of a Term

The effective degrees of freedom of a built term, computed from the
pieces a fitted model supplies. The counting rule follows the term's
penalty. An unpenalized term counts its coefficients exactly. A term
with a smooth penalty uses the trace of \\(H + S)^{-1} H\\ over its
block, where \\H\\ is the term's unpenalized curvature (the weighted
crossproduct of its design block at the fit) and \\S\\ is the penalty's
Hessian in the coefficients at the estimated hyperparameters; without a
penalty the trace reduces to the coefficient count, and as the penalty
grows it falls toward the penalty's null space. A term with a non-smooth
penalty counts its nonzero coefficients, which for the lasso is the
unbiased estimator of its degrees of freedom \[Zou, Hastie & Tibshirani
(2007)\].

## Usage

``` r
edf(term, coef = NULL, hessian = NULL, theta = NULL, tol = 1e-08, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- coef:

  The fitted coefficients of the term's block.

- hessian:

  The unpenalized curvature of the fit restricted to the term's block, a
  \\k \times k\\ matrix; required only for a smooth penalized term.

- theta:

  The estimated hyperparameters of the term's penalty, as a named list;
  required only for a smooth penalized term.

- tol:

  The threshold below which a coefficient counts as zero for a
  non-smooth penalty.

- ...:

  Passed to methods.

## Value

A single number.

## References

Zou, H., Hastie, T. and Tibshirani, R. (2007). On the "degrees of
freedom" of the lasso. *The Annals of Statistics*, 35(5), 2173–2192.

## Examples

``` r
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
built <- term_build(ridge(~ x1 + x2), dd)
H <- crossprod(term_matrix(built))
edf(built, coef = c(0.5, -0.2), hessian = H, theta = list(sigma = 2))
#> [1] 1.976583
```
