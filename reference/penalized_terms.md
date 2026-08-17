# What the Penalized Terms Share

The input handling, the standardization and the prediction of the five
penalized terms, documented once. Each of them –
[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp`](https://statmodels7.github.io/modelterms7/reference/mcp.md) –
has a page of its own carrying its formula, its hyperparameters and
where those may lie. Each takes its block as a one-sided formula or as a
numeric matrix and attaches the corresponding penalties7 object to the
block's coefficients at build time, so the hyperparameters, their bounds
and links, the derivatives and the kink set are the penalty's, never
restated by the term.

## Arguments

- x:

  A one-sided formula or a numeric matrix.

- label:

  A single non-empty string prefixed to the coefficient names.

- standardize:

  A single logical: whether to penalize each coefficient on the scale of
  its own column. See the section below.

- sparse:

  Governs the FORMULA route: whether the block is built as a `dgCMatrix`
  through
  [`sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)
  rather than as a dense model matrix. `NULL`, the default, settles it
  at build from the size of the design; `TRUE` and `FALSE` override it.
  A MATRIX input needs no such argument, being kept in whatever storage
  it arrives in. See the section below.

- ...:

  Not used, and reported: an argument named after another penalty's
  hyperparameter is the mistake this catches.

## Value

An object of class
[`PenalizedTerm`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

A formula input goes through the
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) machinery
with the intercept removed (a penalized block does not penalize an
intercept; the model's intercept lives in the parametric block), and its
blueprint records the terms, the factor levels and the contrasts,
exactly as
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
does. The exception is a formula whose intercept is all it has:
`ridge(~1)` is a block of that one column under the penalty, since
removing it would leave no block. That is the form a subformula on
another term's parameter uses, `gamma ~ lasso(~1)` saying that the
parameter itself carries a lasso, there being no parametric block of its
own to hold an unpenalized intercept. A matrix input is used as given,
and its columns are named after the matrix's own column names, or
numbered when it has none.

Prediction for a matrix input re-evaluates the expression that produced
the matrix in the new data, and ONLY there, so the intended use is a
matrix column of the model data frame (`dd$R <- R`; then `ridge(R)` in
the formula): a subset of the data then carries the matching rows. A
free-standing matrix from the calling environment builds, since its
value was captured, but prediction is rejected – resolving it outside
the new data would silently reuse the build-time rows.

[`term_smooth`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
is `TRUE` for `ridge` and `FALSE` for `lasso`, `enet`, `scad` and `mcp`,
read from each penalty's kink set.

## Sparse storage

A MATRIX input is kept in whatever storage it arrives in, so a
`dgCMatrix` passed to any of the five stays one and nothing has to be
said. `sparse = TRUE` governs the FORMULA route, where the model matrix
would otherwise be built dense whatever the columns look like: it goes
through
[`sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
which BUILDS the block sparse rather than building a dense one and
compressing it.

It pays where the formula carries a factor of many levels, whose
indicator columns hold one non-zero per row – `lasso(~ 0 + g)` over
hundreds of groups is the case. On numeric covariates the block is dense
whatever is asked for, and the sparse storage costs more than it saves.
Left `NULL`, which is the default, the storage is settled at build by
[`.resolve_sparse`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md):
the dense indicator part holds `n` times its column count in cells
against one non-zero per row, and the two routes cross at about \\10^5\\
of those cells.

Standardization does not interfere: it is a diagonal map on the PENALTY
and never an operation on the design, so a sparse block stays sparse
under it.

## Standardization

A hyperparameter is comparable across coordinates only where the
coordinates share a scale: without `standardize` a lasso penalizes a
column measured in metres more than the same column measured in
kilometres, and a reader of \\\lambda\\ has no way to know.

`standardize = TRUE` divides each coefficient by the standard deviation
of its own column, and it does so through the penalty's diagonal map
rather than by touching the design. With \\z_j = x_j/s_j\\ the
coefficient satisfies \\\beta\_{z,j} = s_j\beta\_{x,j}\\, so

\$\$\lambda\sum_j \lvert\beta\_{z,j}\rvert = \lambda\sum_j
s_j\lvert\beta\_{x,j}\rvert = \rho(S\beta_x), \qquad S =
\mathrm{diag}(s),\$\$

which is the standardized penalty read on the original scale. Three
things follow. The design is never rescaled, so a sparse block stays
sparse; \\\lambda\\ stays one number and the coefficients are already on
the scale the data came in, with nothing to map back; and centring,
which is what would destroy sparsity, is not needed, the fit being
invariant to a translation of a penalized column wherever an intercept
is free.

The spread is computed from the built block and frozen in the blueprint,
so the same term standardizes identically in every equation of a
distributional model and does not move with the working weights of a
fit. A constant column takes \\s_j = 1\\.
[`print`](https://rdrr.io/r/base/print.html) shows the values, a number
that changes the meaning of \\\lambda\\ having to be legible.

For SCAD and MCP the diagonal map is not a rescaling of \\\lambda\\
alone: substituting \\s_j\beta_j\\ gives \\\lambda_j = \lambda s_j\\ AND
\\a_j = a/s_j\\ (or \\\gamma_j = \gamma/s_j\\), a composition of both
hyperparameters per coordinate, which the map expresses exactly.

[`random`](https://statmodels7.github.io/modelterms7/reference/random.md)
does not standardize and takes no such argument. Its columns are
grouping indicators and its penalty is a variance component with a
meaning of its own; weighting it by the size of the groups would change
the model rather than its parametrization.

## See also

[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
[`s`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`random`](https://statmodels7.github.io/modelterms7/reference/random.md),
[`term_penalty`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
[`edf`](https://statmodels7.github.io/modelterms7/reference/edf.md)

## Examples

``` r
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
built <- term_build(lasso(~ x1 + x2), dd)
term_coef_names(built)
#> [1] "lasso.x1" "lasso.x2"
term_penalty(built)@params
#> [1] "lambda"
term_smooth(built)
#> [1] FALSE

# the same block penalized on a common scale
dd$x3 <- 1000 * dd$x2
term_penalty(term_build(lasso(~ x1 + x3, standardize = TRUE), dd))@map
#> 2 x 2 diagonal matrix of class "ddiMatrix"
#>         [,1]    [,2]
#> [1,] 1.28964       .
#> [2,]       . 1088.56
```
