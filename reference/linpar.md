# Unpenalized Parametric Term

Creates the specification of an unpenalized parametric block: the design
matrix of a one-sided formula, with the usual
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) conventions
for factors, contrasts, interactions and the intercept.

## Usage

``` r
linpar(formula, label = "", sparse = NULL, contrasts = NULL)
```

## Arguments

- formula:

  A one-sided formula, e.g. `~ x1 + x2`.

- label:

  A character string; when non-empty it is prefixed to the coefficient
  names as `label.name`.

- sparse:

  Whether to build the block as a `dgCMatrix`. `NULL`, the default,
  settles it from the design. See the section below.

- contrasts:

  The contrasts for the formula's factors, as a named list of the kind
  [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)'s
  `contrasts.arg` takes. `NULL`, the default, leaves them to the
  session's `options("contrasts")`.

## Value

An object of class
[`LinparTerm`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The block is the model matrix \\X\\ of the formula and the term's
contribution to the predictor is linear in its coefficients,

\$\$\eta = X\beta,\$\$

with no penalty attached, so all \\p = \operatorname{ncol}(X)\\
coefficients are free and
[`edf`](https://statmodels7.github.io/modelterms7/reference/edf.md)
counts every one of them.

[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
collects the bare covariates of a model formula into one term of this
kind, so `y ~ x1 + x2` and `y ~ linpar(~ x1 + x2)` produce the same
block; the explicit constructor exists for callers who want several
parametric blocks with distinct labels.

Building the term records a blueprint: the terms object, the factor
levels and the contrasts.
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reapplies the mapping through that blueprint, so a factor column in new
data is encoded against the levels seen at build time, and a level the
blueprint does not know is rejected rather than re-encoded. Missing
values are propagated (`na.pass`), never dropped, so the block stays
row-aligned with the response.

## Sparse storage

`sparse = TRUE` builds the block through
[`sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
which BUILDS it sparse rather than building a dense matrix and
compressing it – the second would cost the memory the choice exists to
avoid. Measured at 20000 rows and a factor of 1000 levels, 0.002 s and
1.8 MB against 0.100 s and 161.5 MB, the numbers identical; and a design
that would be 32 GB dense builds in 0.02 s and 19 MB, which is what says
there is no dense intermediate.

It pays where the formula carries a FACTOR OF MANY LEVELS, whose
indicator columns hold one non-zero per row. On numeric covariates the
block is dense whatever is asked for, and the sparse storage then costs
more than it saves. `sparse = NULL`, the default, settles it at build
from the design: the dense indicator part holds `n` times its column
count in cells against one non-zero per row, and the two routes cross at
about \\10^5\\ of those cells, which is the rule
[`.resolve_sparse`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md)
applies. `TRUE` and `FALSE` override it. The storage that was settled is
part of the blueprint, so
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
builds new data the same way.

## See also

[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`enet`](https://statmodels7.github.io/modelterms7/reference/enet.md)

## Examples

``` r
dd <- data.frame(x = 1:4, g = factor(c("a", "a", "b", "b")))
built <- term_build(linpar(~ x + g), dd)
term_matrix(built)
#>   (Intercept) x gb
#> 1           1 1  0
#> 2           1 2  0
#> 3           1 3  1
#> 4           1 4  1
term_coef_names(built)
#> [1] "(Intercept)" "x"           "gb"         
```
