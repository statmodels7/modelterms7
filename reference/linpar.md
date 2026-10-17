# Unpenalized Parametric Term

Specifies an unpenalized parametric block: the model matrix of a
one-sided formula, with the usual
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
conventions for factors, contrasts, interactions and the intercept. It
is the plainest term in the package, and
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
collects a formula's bare covariates into it.

The term's contribution to the predictor is \\X\beta\\, with no penalty,
so every column costs one degree of freedom and
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
returns the column count exactly.

## Usage

``` r
linpar(formula, label = "", sparse = NULL, contrasts = NULL)
```

## Arguments

- formula:

  A one-sided formula, `~ x1 + x2`. A two-sided formula is rejected with
  the error `"'formula' must be one-sided, e.g. ~ x1 + x2."`, and so is
  any object that is not a formula. Its environment is kept and used for
  symbols that the data do not carry. A formula with no columns at all,
  `~ 0`, is rejected at build, when the empty set of column names is
  assigned to the `coef_names` property.

- label:

  A single character string, `""` by default. When non-empty it is
  prefixed to every coefficient name as `label.name`, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) uses it as
  the title. Any other value is rejected.

- sparse:

  `TRUE` to build a `dgCMatrix`, `FALSE` for a base matrix, or `NULL`,
  the default, to settle it from the design by the rule above. Any other
  value signals the error
  `"'sparse' in 'linpar' must be TRUE, FALSE, or NULL to settle it from the design."`.

- contrasts:

  A named list of contrasts for the formula's factors, of the kind that
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)'s
  `contrasts.arg` takes, or `NULL`, the default, for the session's
  `options("contrasts")`. A value that is not a list is rejected.

## Value

An unbuilt
[`LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md):
a specification in which `X`, `coef_names` and `blueprint` are empty
until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills them, and `penalty` is `NULL` permanently.

## The block, and what a build records

Building runs
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) with
`na.action = na.pass` and `drop.unused.levels = FALSE`, so a row with a
missing covariate keeps its place and the block stays aligned with the
response, and a factor level present in the data but used by no row
still gets its column.

The blueprint records the terms object, the factor levels
([`stats::.getXlevels()`](https://rdrr.io/r/stats/checkMFClasses.html)),
the contrasts actually used, the storage chosen at build and the
`drop_intercept` setting.
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reapplies each of them, so a factor column at new rows is encoded
against the levels seen at build time. A level that the blueprint does
not record is rejected by
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) with
`"factor g has new levels zz"`.

## Several parametric blocks

`y ~ x1 + x2` and `y ~ linpar(~ x1 + x2)` give the same columns,
`(Intercept)`, `x1` and `x2`, but the second splits them into two
blocks: the block that
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
builds implicitly, which holds the intercept, and the block of
`linpar()`, which holds `x1` and `x2`. The explicit constructor is
therefore intended for callers who want more than one parametric block
with distinct labels, or who want to set `sparse` or `contrasts` on a
block. Arguments for the block that
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
builds implicitly go through that function's own `linpar` argument.

The intercept belongs to the formula. When the formula carries one, a
`linpar()` written in it drops its own column `(Intercept)` after the
model matrix is built, so its factors keep their contrasts and, together
with the implicit block, `y ~ linpar(~ x + g)` gives the columns of
`lm(y ~ x + g)`. Under `0 +` the block keeps its intercept.

## Sparse storage

`sparse = TRUE` builds through
[`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
which builds the block sparse. Building a dense matrix and compressing
it would use the memory that the sparse route is meant to avoid.

Sparse storage pays where the formula carries a factor of many levels,
whose indicator columns hold one non-zero per row. On numeric covariates
the block has few zero entries, and sparse storage then costs more than
it saves.

`sparse = NULL`, the default, settles the storage at build from the size
of the indicator part: the block is built sparse when
`n * ncol_ind > 1e5`, where `ncol_ind` is an upper bound on the number
of columns coming from factors
([`.indicator_cols()`](https://statmodels7.github.io/modelterms7/reference/dot-indicator_cols.md)).

The settled storage is part of the blueprint, so
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
returns the same kind of block at new rows.

## See also

[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
and
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md)
for the penalized blocks;
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
for the penalized structures;
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
which builds one of these implicitly;
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
and
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md).

## Examples

``` r
dd <- data.frame(x = 1:8, g = factor(rep(c("a", "b", "c", "d"), 2)))

built <- term_build(linpar(~ x + g), dd)
term_matrix(built)
#>   (Intercept) x gb gc gd
#> 1           1 1  0  0  0
#> 2           1 2  1  0  0
#> 3           1 3  0  1  0
#> 4           1 4  0  0  1
#> 5           1 5  0  0  0
#> 6           1 6  1  0  0
#> 7           1 7  0  1  0
#> 8           1 8  0  0  1
term_coef_names(built)
#> [1] "(Intercept)" "x"           "gb"          "gc"          "gd"         

# Unpenalized: every column costs one degree of freedom.
c(npar = term_npar(built), edf = edf(built))
#> npar  edf 
#>    5    5 

# A label prefixes the names, which is how two blocks stay apart.
term_coef_names(term_build(linpar(~ x, label = "lin"), dd))
#> [1] "lin.(Intercept)" "lin.x"          

# Contrasts are recorded at build and reapplied at prediction.
bc <- term_build(linpar(~ g, contrasts = list(g = "contr.sum")), dd)
term_coef_names(bc)
#> [1] "(Intercept)" "g1"          "g2"          "g3"         

# A missing covariate keeps its row, so the block stays aligned.
term_matrix(term_build(linpar(~ x),
                       data.frame(x = c(1, NA, 3))))
#>   (Intercept)  x
#> 1           1  1
#> 2           1 NA
#> 3           1  3

# Storage is settled from the design: a small factor stays dense.
small <- data.frame(g = factor(rep(1:3, 10)))
class(term_matrix(term_build(linpar(~ g), small)))
#> [1] "matrix" "array" 

# Two hundred levels over two thousand rows is 4e5 cells, so sparse.
set.seed(1)
big <- data.frame(g = factor(sample(1:200, 2000, TRUE)))
class(term_matrix(term_build(linpar(~ g), big)))
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"


# The data are simulated from a known truth, so the estimates
# below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(1)
  fd <- data.frame(x = runif(200, -2, 2),
                   g = factor(sample(c("a", "b"), 200, TRUE)))
  fd$y <- 1 + 2 * fd$x + 0.5 * (fd$g == "b") + rnorm(200, sd = 0.4)
  ft <- statmodels7::statmod(y ~ x + g,
                             distributions7::gaussian1_distrib(), fd)
  # The true values are intercept 1, slope 2, and 0.5 for the second level
  round(coef(ft)$mu, 2)
}
#> (Intercept)           x          gb 
#>        1.04        2.00        0.45 
```
