# Unpenalized Parametric Term

Specifies an unpenalized parametric block: the model matrix of a
one-sided formula, with the usual
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
conventions for factors, contrasts, interactions and the intercept. It
is the plainest term in the package, and the one
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
collects a formula's bare covariates into.

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

  A one-sided formula, `~ x1 + x2`. A two-sided formula throws
  `"'formula' must be one-sided, e.g. ~ x1 + x2."`, and anything that is
  not a formula throws. Its environment is kept and used for symbols the
  data do not carry. A formula with no columns at all, `~ 0`, fails in
  the class validator at build time.

- label:

  A single character string, `""` by default. When non-empty it is
  prefixed to every coefficient name as `label.name`, and it is the
  title [`plot()`](https://rdrr.io/r/graphics/plot.default.html) uses.
  Anything that is not one string throws.

- sparse:

  `TRUE` to build a `dgCMatrix`, `FALSE` for a base matrix, or `NULL`,
  the default, to settle it from the design by the rule above. Anything
  else throws
  `"'sparse' in 'linpar' must be TRUE, FALSE, or NULL to settle it from the design."`.

- contrasts:

  A named list of contrasts for the formula's factors, of the kind
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)'s
  `contrasts.arg` takes, or `NULL`, the default, for the session's
  `options("contrasts")`. Anything that is not a list throws.

## Value

An unbuilt
[`LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md):
a specification, with `X`, `coef_names` and `blueprint` empty until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills them, and `penalty` `NULL` permanently.

## The block, and what a build records

Building runs
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) with
`na.action = na.pass` and `drop.unused.levels = FALSE`, so a row with a
missing covariate keeps its place and the block stays aligned with the
response, and a factor level present in the data but used by no row
still gets its column.

The blueprint records the terms object, the factor levels
([`stats::.getXlevels()`](https://rdrr.io/r/stats/checkMFClasses.html)),
the contrasts actually used and the storage settled on.
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reapplies all four, so a factor column at new rows is encoded against
the levels seen at build time. A level the blueprint does not know is
rejected by
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) with
`"factor g has new levels zz"`.

## Several parametric blocks

`y ~ x1 + x2` and `y ~ linpar(~ x1 + x2)` produce the same block, so the
explicit constructor is for callers wanting more than one parametric
block with distinct labels, or wanting to set `sparse` or `contrasts` on
it. Arguments for the block
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
builds implicitly go through that function's own `linpar` argument.

## Sparse storage

`sparse = TRUE` builds through
[`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
which builds the block sparse. Building a dense matrix and compressing
it would cost the memory the choice exists to avoid. Measured at 20000
rows and a factor of 1000 levels, the two routes give identical numbers
at 0.007 s and 1.8 MB against 0.164 s and 161.3 MB; a design of 20000 by
19014 that would be 3.0 GB dense builds sparse in 0.3 s and 3.1 MB.

It pays where the formula carries a **factor of many levels**, whose
indicator columns hold one non-zero per row. On numeric covariates the
block is dense whatever is asked for, and the sparse storage then costs
more than it saves.

`sparse = NULL`, the default, settles it at build from the size of the
indicator part: `n * ncol_ind > 1e5`, where `ncol_ind` counts the
columns coming from factors
([`.indicator_cols()`](https://statmodels7.github.io/modelterms7/reference/dot-indicator_cols.md)).
Measured over fifteen combinations, building the block and forming its
crossproduct, the two routes cross between \\10^5\\ and \\3 \times
10^5\\ cells and the sparse route then wins by orders: at 20000 rows it
is 1.4 times faster at 15 levels, 14 times at 60 and 445 times at 400.

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
```
