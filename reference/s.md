# Penalized Smooth of One Covariate

A smooth function of one covariate, expanded in a basis7 basis and
penalized for roughness. The default is a cubic B-spline basis under the
Demmler-Reinsch reparametrization, which separates the linear effect
from the nonlinear deviation and turns the roughness penalty into the
identity on the deviation.

As the smoothing parameter grows the fit approaches a straight line, and
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
falls to exactly one.

## Usage

``` r
s(
  x,
  by = NULL,
  k = 10,
  degree = 3,
  basis = NULL,
  linear = TRUE,
  label = NULL,
  lambda = NULL,
  id = NULL,
  sparse = NULL
)
```

## Arguments

- x:

  The covariate, an expression evaluated in the data.

- by:

  An optional factor or numeric variable, given as a bare expression;
  `NULL` by default. See the section above.

- k:

  The basis dimension before reparametrization, `10` by default. It must
  exceed `degree`: a cubic spline needs at least four basis functions,
  and anything smaller throws. The block has `k - 1` columns with
  `linear = TRUE` and `k - 2` without it.

- degree:

  The spline degree, `3` by default, a cubic spline.

- basis:

  An optional basis7 basis used in place of the default B-spline. Its
  range is taken as given, so a basis built on one interval is not
  re-placed on the data's.

- linear:

  Whether the linear effect is carried in the block and left unpenalized
  there, `TRUE` by default.

- label:

  A single non-empty string prefixed to the coefficient names. `NULL`,
  the default, builds one from the covariate: `s(x)`.

- lambda:

  The smoothing parameter, held at the value given and **estimated**
  when left `NULL`, which is the default.

- id:

  A label sharing this smooth's smoothing parameter with those of other
  terms carrying the same one: they are then estimated at a single
  value, so several curves are smoothed together. It is what `id` does
  in mgcv, and what it means best between smooths of the same basis and
  dimension. `NULL`, the default, shares nothing. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md).

- sparse:

  `TRUE`, `FALSE`, or `NULL` to settle it at build. Only a factor `by`
  admits `TRUE`; without one it is refused rather than ignored. See the
  section above.

## Value

An unbuilt
[`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md):
a specification, with `X`, `coef_names`, `blueprint` and `penalty` empty
until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills them.

## The block and its penalty

The block has one column for the linear effect, centered and scaled,
followed by the reparametrized basis, so `s(x, k = 8)` gives seven
columns named `s(x).lin`, `s(x).z1` ... `s(x).z6`. That ordering is what
the penalty reads: it is the quadratic penalty of \\\mathrm{diag}(0, 1,
\dots, 1)\\, rank deficient by exactly one, so the linear effect is
unpenalized and the deviation is shrunk toward zero.

Two consequences a reader of a fit needs. At a large smoothing parameter
the fit tends to a straight line, so
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
runs from `k - 1` down to 1 and never to 0. And the linear column is
orthogonal to the rest over the observed covariate, so the linear and
the nonlinear parts of a fitted smooth are separately readable.

`linear = FALSE` drops that first column, and the penalty is then the
identity over the whole block, of full rank.

## The construction is empirical

The Demmler-Reinsch transform (Demmler and Reinsch, 1975; used for
effect selection by Bach and Klein, 2024) takes the inner product **at
the observed covariate values**, so it is computed when the term is
built and stored in the blueprint. Prediction is the parent basis
evaluated at the new points and multiplied by that same transform, so
the separation of the linear from the nonlinear part holds at new rows
as it does at old.

Rebuilding on other rows instead would place the knots on their range
and compute another transform. Measured on 80 points, predicting on the
first ten agrees with those rows of the original block exactly and
rebuilding differs by 2.85.

## Varying the smooth by another variable

A **factor** `by` gives one smooth per level: the block is the smooth
multiplied by each level's indicator, and the penalty is the same matrix
repeated blockwise, so one smoothing parameter governs every level.
`s(x, k = 5, by = g)` over a four-level factor has 16 columns.

A **numeric** `by` gives a varying-coefficient term: the smooth
multiplies that variable, and the fitted function is the coefficient of
`by` as it changes with the covariate.

## Sparse storage

A factor `by` is the one place a smooth's block can be sparse, each row
sitting in the block of its own level and nowhere else, a density of
\\1/m\\. `sparse = TRUE` builds it that way instead of building the
dense matrix and compressing it: measured at 2000 rows, \\k = 10\\ and
200 levels, 0.35 MB against 28.93 MB with the numbers identical.

`sparse = NULL`, the default, settles it at build from the size of the
block through
[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md),
the dense form holding \\n m k\\ cells against \\n k\\ non-zeros.

An explicit `TRUE` is **refused** without a factor `by`, and the message
says why: the basis is dense by construction and a numeric `by` merely
multiplies it, so there would be nothing to build on.

The block alone is sparse. The penalty of a factor `by` is the same
matrix repeated blockwise and penalties7 returns it dense, 25.92 MB at
those sizes; that is a property of that package's contract.

## References

Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
smoothing. *Numerische Mathematik*, 24, 375–382.

Bach, P. and Klein, N. (2024). Bayesian effect selection in additive
models with an application to time-to-event data.

## See also

[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) for
several covariates,
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
for a grouped effect,
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) for
a parametric nonlinear shape,
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what a fitted smooth spends.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(80)), g = factor(rep(letters[1:4], 20)))
dd$y <- sin(2 * pi * dd$x) + rnorm(80, sd = 0.2)

# k = 8 gives seven columns: the linear effect and six deviations.
b <- term_build(s(x, k = 8), dd)
term_coef_names(b)
#> [1] "s(x).lin" "s(x).z1"  "s(x).z2"  "s(x).z3"  "s(x).z4"  "s(x).z5"  "s(x).z6" 

# The penalty is diag(0, 1, ..., 1): the linear column is free.
penalties7::penalty_matrix(term_penalty(b), list(lambda = 1))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]    0    0    0    0    0    0    0
#> [2,]    0    1    0    0    0    0    0
#> [3,]    0    0    1    0    0    0    0
#> [4,]    0    0    0    1    0    0    0
#> [5,]    0    0    0    0    1    0    0
#> [6,]    0    0    0    0    0    1    0
#> [7,]    0    0    0    0    0    0    1

# The linear column really is the linear effect, and is orthogonal to
# the rest over the observed covariate.
X <- term_matrix(b)
cor(X[, 1], dd$x)
#> [1] 1
max(abs(crossprod(X[, 1], X[, -1])))
#> [1] 4.418688e-14

# So edf runs from k - 1 down to one, not to zero.
H <- crossprod(X)
cf <- rnorm(ncol(X))
vapply(c(1e-8, 1, 1e12),
       function(l) edf(b, coef = cf, hessian = H, theta = list(lambda = l)),
       numeric(1))
#> [1] 7.000000 4.191177 1.000000

# A factor `by` is one smooth per level under one smoothing parameter.
bf <- term_build(s(x, k = 5, by = g), dd)
c(npar = term_npar(bf), levels = nlevels(dd$g))
#>   npar levels 
#>     16      4 

# The transform is computed on the data and reapplied, never rebuilt.
max(abs(term_predict(b, dd[1:10, ]) - X[1:10, ]))
#> [1] 0
max(abs(term_matrix(term_build(s(x, k = 8), dd[1:10, ])) - X[1:10, ]))
#> [1] 3.628495

# Sparsity needs a factor `by`, and says so when there is none.
try(term_build(s(x, k = 5, sparse = TRUE), dd))
#> Error : 'sparse' has nothing to build on here: a smooth's basis is dense by
#>   construction. Sparsity comes from a FACTOR 'by', whose indicators put each
#>   row in the block of its own level.


# Fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(4)
  fd <- data.frame(z = sort(runif(200, -3, 3)))
  fd$y <- sin(fd$z) + rnorm(200, sd = 0.3)
  ft <- statmodels7::statmod(y ~ s(z),
                             distributions7::gaussian1_distrib(), fd)
  # the smoothing parameter is chosen by REML, and the fit follows sin()
  round(c(edf = sum(ft@edf$edf),
          rmse = sqrt(mean((fitted(ft) - sin(fd$z))^2))), 3)
}
#>   edf  rmse 
#> 8.510 0.027 
```
