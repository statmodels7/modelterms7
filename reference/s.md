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
  smoother = basis7::bspline_smooth(),
  by = NULL,
  by_hyper = c("shared", "level"),
  hyper = NULL,
  id = NULL,
  label = NULL,
  sparse = NULL,
  ...
)
```

## Arguments

- x:

  The covariate, an expression evaluated in the data.

- smoother:

  How the smooth is built: a basis7 smoother, which carries the basis,
  the roughness penalty, the null space and the reparametrization
  together.
  [`basis7::bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.html)
  is the default and is the construction this term has always used;
  [`basis7::fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.html)
  is the periodic one and
  [`basis7::legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.html)
  the global polynomial one.

- by:

  An optional factor or numeric variable, given as a bare expression;
  `NULL` by default. See the section above.

- by_hyper:

  With a factor `by`, whether the levels share one smoothing parameter
  (`"shared"`, the default) or carry one each (`"level"`). See the
  section above. Rejected without a factor `by`, which has no levels to
  give one each.

- hyper:

  The hyperparameters of the smoother's penalty to hold, as a named
  numeric vector such as `c(lambda = 2)`. What names there are depends
  on the penalty; anything left out is **estimated**, which is the
  default. Under `by_hyper = "level"` a name is qualified by the level,
  `c(lambda.a = 2)`.

- id:

  A label sharing this smooth's smoothing parameter with those of other
  terms carrying the same one: they are then estimated at a single
  value, so several curves are smoothed together. It is what `id` does
  in mgcv, and what it means best between smooths of the same basis and
  dimension. `NULL`, the default, shares nothing. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md).

- label:

  A single non-empty string prefixed to the coefficient names. `NULL`,
  the default, builds one from the covariate: `s(x)`.

- sparse:

  `TRUE`, `FALSE`, or `NULL` to settle it at build. Only a factor `by`
  admits `TRUE`; without one it is refused rather than ignored. See the
  section above.

- ...:

  Not used, and accepted only so that an argument `s()` no longer takes
  is reported with its replacement rather than as R's own "unused
  argument", which names the argument and not what to write.

## Value

An unbuilt
[`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md):
a specification, with `X`, `coef_names`, `blueprint` and `penalty` empty
until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills them.

## The block and its penalty

The block has one column for the linear effect, centered and scaled,
followed by the reparametrized basis, so `s(x, bspline_smooth(k = 8))`
gives seven columns named `s(x).lin`, `s(x).z1` ... `s(x).z6`. That
ordering is what the penalty reads: it is the quadratic penalty of
\\\mathrm{diag}(0, 1, \dots, 1)\\, rank deficient by exactly one, so the
linear effect is unpenalized and the deviation is shrunk toward zero.

Two consequences a reader of a fit needs. At a large smoothing parameter
the fit tends to a straight line, so
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
runs from `k - 1` down to 1 and never to 0. And the linear column is
orthogonal to the rest over the observed covariate, so the linear and
the nonlinear parts of a fitted smooth are separately readable.

`bspline_smooth(null_space = "drop")` drops that first column, and the
penalty is then the identity over the whole block, of full rank.

## A penalty of your own

The smoother's `penalty` argument replaces the roughness matrix with a
penalties7 penalty of its own, built at the coefficient count only the
data settle:
`s(x, bspline_smooth(k = 20, penalty = penalties7::lasso_penalty))`. A
penalties7 constructor passes bare, and anything else is a function of
the count.

What it buys is that the Demmler-Reinsch coordinates are ordered from
the smoothest to the most wiggly and the roughness penalty on them is
the identity, so an \\\ell_1\\ penalty takes whole directions of
wiggliness to exactly zero and chooses the smooth's effective dimension.
Measured at \\n = 300\\ with \\k = 20\\, the coordinates surviving as
the smoothing parameter grows are 15, 10, 4, 2, 1 and 0 of 18, and the
fit at two of them is as close to the truth as the fit at all eighteen.
A heavy-tailed penalty is the robust reading of the same block, and a
structured one estimates the correlation of the coefficients rather than
fixing it.

Two things follow, and neither is a detail. The penalty covers **the
penalized coordinates alone**: a roughness matrix carries a zero row for
the linear column and leaves it free by its own arithmetic, where a
separable penalty has no such row and would shrink it, so the free
columns are outside the entry it declares. And a penalty **with a kink**
has no second derivative at zero, which a marginal criterion needs, so
its smoothing parameter is chosen by a path over `sparse_criterion` –
BIC by default – and not by REML, whatever `outer_criterion` says. That
is a change in what the number means and not only in how it is computed.

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
repeated blockwise. `s(x, bspline_smooth(k = 5), by = g)` over a
four-level factor has 16 columns.

`by_hyper` says whether the levels are smoothed **together or apart**.
`"shared"`, the default, estimates one smoothing parameter for all of
them, which is what this package has always done. `"level"` declares one
penalty per level through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
so each is estimated on its own; it is what mgcv does by default, and
what `id` there undoes. Which to want is a question about the data
rather than about the code: levels that differ in how wiggly they are,
or in how much of them there is, ask for different amounts of smoothing.
Measured on three levels of one curve at amplitudes 1, 0.15 and 2.5, the
separately estimated smoothing parameters are 2.03, 111.9 and 0.468 – a
factor of 239 apart – and the effective degrees of freedom fall from
25.5 to 21.7 as the flattest level is smoothed away.

Under `"level"` a held hyperparameter is qualified by the level it
belongs to, `hyper = c(lambda.a = 2)`, and so is an `id`.

A **numeric** `by` gives a varying-coefficient term: the smooth
multiplies that variable, and the fitted function is the coefficient of
`by` as it changes with the covariate. It has no levels, so
`by_hyper = "level"` is rejected rather than read as `"shared"`.

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
b <- term_build(s(x, basis7::bspline_smooth(k = 8)), dd)
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
#> [1] 8.992806e-14

# So edf runs from k - 1 down to one, not to zero.
H <- crossprod(X)
cf <- rnorm(ncol(X))
vapply(c(1e-8, 1, 1e12),
       function(l) edf(b, coef = cf, hessian = H, theta = list(lambda = l)),
       numeric(1))
#> [1] 7.000000 4.191177 1.000000

# A factor `by` is one smooth per level under one smoothing parameter.
bf <- term_build(s(x, basis7::bspline_smooth(k = 5), by = g), dd)
c(npar = term_npar(bf), levels = nlevels(dd$g))
#>   npar levels 
#>     16      4 
length(term_penalties(bf))
#> [1] 1

# by_hyper = "level" gives them one each: one entry per level, over that
# level's own columns.
bl <- term_build(s(x, basis7::bspline_smooth(k = 5), by = g,
                   by_hyper = "level"), dd)
vapply(term_penalties(bl), function(e) e$name, character(1))
#> [1] "a" "b" "c" "d"
vapply(term_penalties(bl), function(e) range(e$index), integer(2))
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    5    9   13
#> [2,]    4    8   12   16

# A penalty factory covers the PENALIZED coordinates: the free linear
# column is left out, a separable penalty having no zero row with which
# to leave it alone.
bp <- term_build(s(x, basis7::bspline_smooth(
  k = 5, penalty = penalties7::lasso_penalty)), dd)
term_penalties(bp)[[1]]$index
#> [1] 2 3 4
term_penalties(bp)[[1]]$penalty@params
#> [1] "lambda"

# The transform is computed on the data and reapplied, never rebuilt.
max(abs(term_predict(b, dd[1:10, ]) - X[1:10, ]))
#> [1] 0
max(abs(term_matrix(term_build(s(x, basis7::bspline_smooth(k = 8)), dd[1:10, ])) - X[1:10, ]))
#> [1] 2.849289

# Sparsity needs a factor `by`, and says so when there is none.
try(term_build(s(x, basis7::bspline_smooth(k = 5), sparse = TRUE), dd))
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
