# Effective Degrees of Freedom of a Term

The effective degrees of freedom of a built term, computed from the
pieces a fitted model supplies. The counting rule follows the penalties
the term declares through
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
and applies to each of them over the parameters it covers.

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
  \\k \times k\\ matrix; required whenever some parameter is not under a
  kinked penalty.

- theta:

  The estimated hyperparameters. For a term carrying one penalty, that
  penalty's hyperparameters as a named list; for a term carrying
  several, a list of such lists keyed by the penalty names
  [`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
  gives.

- tol:

  The threshold below which a coefficient counts as zero for a
  non-smooth penalty.

- ...:

  Passed to methods.

## Value

A single number.

## Details

A parameter no penalty reaches counts one, exactly. A parameter under a
**non-smooth** penalty counts one when it is away from zero and nothing
when it is at it, which for the lasso is the unbiased estimator of its
degrees of freedom \[Zou, Hastie & Tibshirani (2007)\]. The remaining
parameters – those unpenalized and those under a **smooth** penalty –
are counted together by the trace of \\(H + S)^{-1} H\\ over the
sub-block they occupy, where \\H\\ is the term's unpenalized curvature
there (the weighted crossproduct of its design block at the fit) and
\\S\\ carries each smooth penalty's Hessian in the coefficients at the
estimated hyperparameters, placed at the parameters that penalty covers
and zero elsewhere. An unpenalized parameter contributes a zero row and
column to \\S\\, so the trace returns its one; as a penalty grows the
trace falls toward the dimension of its null space.

The rules compose because they partition the term's parameters, and each
reduces to what the term reported before when one penalty covers the
whole block: the trace over every column for a smooth penalty, the
nonzero count for a kinked one, the coefficient count for none.

`hessian` is asked for over the whole block, and is used at the rows and
columns the trace runs over. It is not needed at all when every penalty
is kinked, since the count is then read from `coef` alone.

## References

Zou, H., Hastie, T. and Tibshirani, R. (2007). On the "degrees of
freedom" of the lasso. *The Annals of Statistics*, 35(5), 2173–2192.

## See also

[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
[`term_penalty`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
[`term_smooth`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)

## Examples

``` r
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
built <- term_build(ridge(~ x1 + x2), dd)
H <- crossprod(term_matrix(built))
edf(built, coef = c(0.5, -0.2), hessian = H, theta = list(lambda = 0.25))
#> [1] 1.976583
```
