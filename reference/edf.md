# Effective Degrees of Freedom of a Term

Counts what a built term spends, given the coefficients a fit reached,
the curvature at that fit and the hyperparameters that were estimated.
An unpenalized parameter costs one; a parameter under a penalty costs
less, and how much less depends on the kind of penalty. The three rules
partition the term's parameters, so a term carrying several penalties is
counted piece by piece and the pieces add.

The one method is registered on
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
and reads the penalties the term declares through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
so a term class written outside the package is counted by the same rule
with nothing to register.

## Usage

``` r
edf(term, coef = NULL, hessian = NULL, theta = NULL, tol = 1e-08, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).
  An unbuilt one throws.

- coef:

  The fitted coefficients of the term's block, a numeric vector of
  length `term_npar(term)`. Any other length throws with the required
  length named. `NULL` is allowed for an unpenalized term.

- hessian:

  The unpenalized curvature of the fit restricted to the term's own
  block, a \\k \times k\\ matrix with \\k\\ = `term_npar(term)`,
  typically the weighted crossproduct of the block at the fitted
  weights. Any other shape throws. `NULL` is allowed unless the trace
  runs.

- theta:

  The estimated hyperparameters. A term carrying one penalty takes that
  penalty's own named list, `list(lambda = 0.25)`. A term carrying
  several takes a list of such lists, keyed by the names
  [`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
  gives; an unkeyed list throws, naming them. The two spellings are told
  apart by the value being a list, a hyperparameter never being one.

- tol:

  The magnitude below which a coefficient counts as zero under a
  non-smooth penalty, `1e-8` by default. A proximal step returns exact
  zeros, so the threshold matters only for coefficients a different
  route left small.

- ...:

  Passed to methods.

## Value

A single number, between 0 and `term_npar(term)`. Not an integer in
general: only the unpenalized and the non-smooth rules give whole
numbers.

## The three rules

**A parameter no penalty reaches counts one**, exactly. A term with no
penalties at all therefore returns `term_npar(term)` and reads none of
the other arguments.

**A parameter under a non-smooth penalty counts one when it is away from
zero and nothing when it is at it.** For the lasso that count is an
unbiased estimator of the degrees of freedom (Zou, Hastie and
Tibshirani, 2007). A penalty is treated as non-smooth when
[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
reports any point at a probe value of its hyperparameters, so lasso,
SCAD, MCP and the elastic net go here and ridge does not.

**Everything else is counted together by one trace.** Over the
parameters left, meaning those unpenalized and those under a smooth
penalty,

\$\$\mathrm{edf} = \mathrm{tr}\\(H + S)^{-1} H\\,\$\$

where \\H\\ is `hessian` restricted to those rows and columns and \\S\\
carries each smooth penalty's Hessian in the coefficients, evaluated at
`coef` and at the estimated hyperparameters and placed at the parameters
that penalty covers. An unpenalized parameter contributes a zero row and
column to \\S\\, so the trace returns its one. As a smoothing parameter
grows the trace falls toward the dimension of the penalty's null space:
for [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md)
that limit is 1, the straight line the Demmler-Reinsch penalty leaves
free.

Each rule reduces to what a term reported before it carried more than
one penalty: the trace over every column for a single smooth penalty,
the nonzero count for a single kinked one, the coefficient count for
none.

## Which arguments are needed when

`hessian` is asked for over the whole block and read at the rows and
columns the trace runs over. It is never needed when every parameter
sits under a kinked penalty, the count being read from `coef` alone, and
`coef` is never needed when the term carries no penalty. Everything the
rule in force does not reach may be left `NULL`; leaving out something
it does reach throws with the missing arguments named.

## References

Zou, H., Hastie, T. and Tibshirani, R. (2007). On the "degrees of
freedom" of the lasso. *The Annals of Statistics*, 35(5), 2173–2192.

## See also

[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the entries the count runs over,
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
for the ceiling,
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
for whether a term's penalized objective is twice differentiable, and
[`penalties7::penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.html)
for the \\S\\ block.

## Examples

``` r
set.seed(2)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))

# No penalty: the coefficient count, and nothing else is read.
edf(term_build(linpar(~ x1 + x2), dd))
#> [1] 3

# A ridge spends between 2 and 0 as lambda grows.
b <- term_build(ridge(~ x1 + x2), dd)
H <- crossprod(term_matrix(b))
vapply(c(1e-6, 0.25, 10, 1e6),
       function(l) edf(b, coef = c(0.5, -0.2), hessian = H,
                       theta = list(lambda = l)), numeric(1))
#> [1] 2.000000e+00 1.979321e+00 1.412199e+00 4.852976e-05

# And that is the trace, computed apart.
all.equal(edf(b, coef = c(0.5, -0.2), hessian = H, theta = list(lambda = 0.25)),
          sum(diag(solve(H + 0.25 * diag(2), H))))
#> [1] TRUE

# A lasso counts survivors, and needs no hessian.
bl <- term_build(lasso(~ x1 + x2 + x3), dd)
edf(bl, coef = c(1, 0, -2), theta = list(lambda = 1))
#> [1] 2

# A smooth falls from k toward the dimension of its null space, which
# for s() is the one straight line the penalty leaves free.
d2 <- data.frame(x = seq(0, 1, length.out = 60))
bs <- term_build(s(x, k = 8), d2)
Hs <- crossprod(term_matrix(bs))
cf <- rnorm(term_npar(bs))
c(k = term_npar(bs),
  vapply(c(1e-8, 1, 1e10),
         function(l) edf(bs, coef = cf, hessian = Hs,
                         theta = list(lambda = l)), numeric(1)))
#>        k                            
#> 7.000000 7.000000 3.969875 1.000000 

# Two penalties on one term: theta is keyed by the entry names.
d3 <- data.frame(x = runif(40, 0, 5), g = factor(rep(1:4, 10)))
nb <- term_build(nl(~ a * exp(-r * x), a ~ 0 + ridge(~ g),
                    r ~ 0 + lasso(~ g)), d3)
vapply(term_penalties(nb), function(e) e$name, character(1))
#> [1] "a::ridge(~g)" "r::lasso(~g)"
Hn <- crossprod(as.matrix(term_matrix(nb)))
edf(nb, coef = c(1, 2, 0, 3, 0.5, 0, 0.2, 0), hessian = Hn,
    theta = list(`a::ridge(~g)` = list(lambda = 2),
                 `r::lasso(~g)` = list(lambda = 1)))
#> [1] 5.333333
```
