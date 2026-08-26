# Structural Checks on a Model Term

Builds `term` against `data` and runs six checks on the result: that the
build succeeds and returns a two-dimensional numeric block with one row
per observation, that the coefficient names are unique and as numerous
as the columns, that
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
agrees with that count, that
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
answers with a single non-missing logical, that
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
on the same data reproduces the block, and that
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
on a subset of rows returns the corresponding rows of it. One row of the
result per check, printed as it goes and returned invisibly.

The last check is the one worth running. A term is supposed to record
its encoding at build time and reapply it; a term that re-derives the
encoding from whatever rows it is handed passes every other check and
fails this one.

## Usage

``` r
check_term(term, data, verbose = TRUE)
```

## Arguments

- term:

  A term specification: any object inheriting from
  [`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md),
  built or unbuilt. Anything that is not a
  [`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
  throws `"'term' must inherit from 'model_term'."`, and a
  [`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
  is rejected by name. An already-built term is accepted and is rebuilt
  against `data`.

- data:

  A data frame carrying every variable the term names. Anything else
  throws `"'data' must be a data frame."`. Any number of rows is
  accepted; with one row the subset check compares the block against
  itself and is uninformative.

- verbose:

  Print one line per check as it completes, `TRUE` by default, in the
  form ` [OK] subset -- 4 rows`. The rows are returned either way, so
  `verbose = FALSE` is the form to use inside a test.

## Value

Invisibly, a data frame of `check`, `status` and `info`, all character,
one row per check:

- `check`:

  `"build"`, `"names"`, `"npar"`, `"smooth"`, `"reproduce"` and
  `"subset"`, in that order.

- `status`:

  `"OK"` or `"FAILED"`.

- `info`:

  The block's dimensions for `build`, with `", sparse"` appended when it
  is an S4 matrix; `"smooth"` or `"non-smooth"` for `smooth`; the number
  of rows kept for `subset`; the condition message where a check threw;
  and `""` otherwise.

A build that throws gives a single row, `check = "build"` and
`status = "FAILED"`, with the message in `info`; no later check is
attempted.

## The two identities

Write \\X\\ for `term_matrix(term_build(term, data))` and \\S\\ for a
subset of the row indices. The last two checks are

\$\$\mathrm{predict}(\mathrm{term}, \mathrm{data}) = X, \qquad
\mathrm{predict}(\mathrm{term}, \mathrm{data}\[S, \]) = X\[S, \],\$\$

both compared by [`all.equal()`](https://rdrr.io/r/base/all.equal.html)
at a relative tolerance of `1e-12`. The second does not follow from the
first: a term that rebuilds its encoding satisfies the first, the rows
being the same ones, and fails the second as soon as \\S\\ omits a
factor level or narrows the range a basis is placed on.

## Which rows the subset takes

The subset is chosen to make that failure reachable. If any column of
`data` is a factor with two or more levels, every row of its last level
is dropped, so a rebuilt encoding is one column short. Failing that, the
first `nrow(data) %/% 2` rows are taken, which still narrows the range
of a numeric covariate.

[`droplevels()`](https://rdrr.io/r/base/droplevels.html) is applied to
the subset before it is passed on, because that is how new data really
reach
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md):
a factor column there carries only the levels its own rows use, and the
rest are known to the blueprint alone. Without the call a plain row
subset carries the original level set along with it and the check cannot
fail.

## What is not checked

The battery covers the design block. It reads neither the penalty a
penalized term attaches nor the hyperparameters it declares, so
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
and
[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)
are never called.

`check_term()` applies to an additive term and rejects a structural one
by name. A structural term
([`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md),
and the marginal break-point terms) contributes no design columns and
registers no
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
method, so there is nothing for the six checks to read; the message says
so and names what such a term supplies instead, which is
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md),
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
and one of
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
or
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md).
No validator covers those.

## See also

[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
and
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md),
the two generics the checks compare;
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
for the block itself;
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
for reading a whole formula into terms.

## Examples

``` r
dd <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))

# A parametric block passes all six.
res <- check_term(linpar(~ x + g), dd)
#>   [OK] build -- 6 x 4 block
#>   [OK] names
#>   [OK] npar
#>   [OK] smooth -- smooth
#>   [OK] reproduce
#>   [OK] subset -- 4 rows
all(res$status == "OK")
#> [1] TRUE

# What the subset check is testing. The chosen rows drop level "c", and
# the reapplied block keeps the four columns the blueprint recorded.
b  <- term_build(linpar(~ x + g), dd)
nd <- droplevels(dd[c(1, 2, 4, 5), ])
levels(nd$g)
#> [1] "a" "b"
dim(term_predict(b, nd))              # 4 x 4, the blueprint's levels
#> [1] 4 4
dim(model.matrix(~ x + g, nd))        # 4 x 3, what a rebuild gives
#> [1] 4 3

# The same failure for a basis: rebuilding places the knots on the
# narrower range, so the columns are different functions of x.
d2  <- data.frame(x = seq(0, 1, length.out = 40))
bs  <- term_build(s(x, k = 6), d2)
X   <- term_matrix(bs)
sub <- 1:20
max(abs(term_predict(bs, d2[sub, , drop = FALSE]) - X[sub, ]))
#> [1] 0
max(abs(term_matrix(term_build(s(x, k = 6), d2[sub, , drop = FALSE])) -
        X[sub, ]))
#> [1] 2.777484

# A sparse block is accepted, and the info column says it is sparse.
check_term(random(~ 1 | g), dd)
#>   [OK] build -- 6 x 3 block, sparse
#>   [OK] names
#>   [OK] npar
#>   [OK] smooth -- smooth
#>   [OK] reproduce
#>   [OK] subset -- 4 rows

# A build that throws gives one row and stops there.
print(check_term(linpar(~ x + nonexistent), dd, verbose = FALSE))
#>   check status                           info
#> 1 build FAILED object 'nonexistent' not found
```
