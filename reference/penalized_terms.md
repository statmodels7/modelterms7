# What the Penalized Terms Share

The input handling, the storage, the standardization and the prediction
of the five penalized terms, documented once.
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
each have a page of their own carrying the penalty's formula, its
hyperparameters and where those may lie.

Each takes its block as a one-sided formula or as a numeric matrix, and
attaches the corresponding penalties7 object to the block's coefficients
at build time. The hyperparameters, their bounds and links, the
derivatives and the kink set are the penalty's and are never restated by
the term.

## Arguments

- x:

  A one-sided formula, such as `~ x1 + x2` or `~ 0 + g`, or a numeric
  matrix, ideally a matrix column of the model data frame. A two-sided
  formula throws.

- label:

  A single non-empty character string prefixed to the coefficient names
  as `label.name`. Each constructor defaults it to its own name, so a
  ridge over `x1` reads `ridge.x1`, and two penalized terms in one
  formula stay apart by their labels.

- standardize:

  A single logical, `FALSE` by default: whether to penalize each
  coefficient on the scale of its own column. See the section above.

- sparse:

  Governs the formula route: `TRUE` builds a `dgCMatrix` through
  [`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
  `FALSE` a dense model matrix, and `NULL`, the default, settles it at
  build from the size of the design. A matrix input needs no such
  argument. See the section above.

- ...:

  Not used, and reported. An argument named after another penalty's
  hyperparameter is the mistake this catches: `scad(~ x, gamma = 1)`
  throws `"'scad' has no argument 'gamma'."` and lists the ones it does
  have.

## Value

An unbuilt
[`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md):
a specification, with `X`, `coef_names`, `blueprint` and `penalty` all
empty until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills them.

## The two inputs

A **formula** goes through the
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
machinery with the intercept removed: a penalized block does not
penalize an intercept, and the model's intercept lives in the parametric
block. Its blueprint records the terms, the factor levels and the
contrasts, exactly as
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
does, so `~ g` and `~ 0 + g` give the same four columns for a four-level
factor.

The exception is a formula whose intercept is all it has. `ridge(~ 1)`
is a block of that one column under the penalty, removing it leaving no
block at all. That is the form a subformula on another term's parameter
uses: `gamma ~ 0 + lasso(~ 1)` says the parameter itself carries a
lasso, there being no parametric block of its own to hold an unpenalized
intercept.

A **matrix** is used as given, and its columns are named after its own
column names, or numbered `1`, `2`, ... when it has none.

## Predicting a matrix input

Prediction re-evaluates the expression that produced the matrix **in the
new data, and only there**. So the intended use is a matrix column of
the model data frame:

    dd$R <- R
    ridge(R)      # in the formula

A subset of `dd` then carries the matching rows of `R`. A free-standing
matrix from the calling environment builds, its value having been
captured, and prediction is refused with a message naming the
expression: resolving it outside the new data would silently reuse the
build-time rows whenever the counts happened to agree.

## Sparse storage

A **matrix** input is kept in whatever storage it arrives in, so a
`dgCMatrix` passed to any of the five stays one and nothing needs to be
said.

`sparse` governs the **formula** route, where the model matrix would
otherwise be built dense whatever the columns look like. `TRUE` goes
through
[`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
which builds the block sparse; building a dense one and compressing it
would cost the memory the choice exists to avoid.

It pays where the formula carries a factor of many levels, whose
indicator columns hold one non-zero per row, and `lasso(~ 0 + g)` over
hundreds of groups is the case. On numeric covariates the block is dense
whatever is asked for, and the sparse storage then costs more than it
saves. Left `NULL`, the default, the storage is settled at build by
[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md):
the dense indicator part holds `n` times its column count in cells
against one non-zero per row, and the two routes cross at about \\10^5\\
of those cells.

Standardization does not interfere. It is a diagonal map on the
**penalty** and never an operation on the design, so a sparse block
stays sparse under it.

## Standardization

A hyperparameter is comparable across coordinates only where the
coordinates share a scale. Without `standardize` a lasso penalizes a
column measured in meters more than the same column measured in
kilometers, and a reader of \\\lambda\\ has no way to tell.

`standardize = TRUE` divides each coefficient by the standard deviation
of its own column, through the penalty's diagonal map. With \\z_j =
x_j/s_j\\ the coefficients satisfy \\\beta\_{z,j} = s_j\beta\_{x,j}\\,
so

\$\$\lambda\sum_j \lvert\beta\_{z,j}\rvert = \lambda\sum_j
s_j\lvert\beta\_{x,j}\rvert = \rho(S\beta_x), \qquad S =
\mathrm{diag}(s),\$\$

which is the standardized penalty read on the original scale. Three
things follow. The design is never rescaled, so a sparse block stays
sparse. \\\lambda\\ stays one number, and the coefficients are already
on the scale the data came in, with nothing to map back. And centering,
which is what would destroy sparsity, is not needed: the fit is
invariant to a translation of a penalized column wherever an intercept
is free.

The spread is computed from the built block by
[`.block_sd()`](https://statmodels7.github.io/modelterms7/reference/dot-block_sd.md)
and frozen in `blueprint$standardize`, so the same term standardizes
identically in every equation of a distributional model and does not
move with the working weights of a fit. A constant column takes \\s_j =
1\\. [`print()`](https://rdrr.io/r/base/print.html) shows the values, a
number that changes the meaning of \\\lambda\\ having to be legible.

For SCAD and MCP the map is not a rescaling of \\\lambda\\ alone.
Substituting \\s_j\beta_j\\ gives \\\lambda_j = \lambda s_j\\ **and**
\\a_j = a/s_j\\ (or \\\gamma_j = \gamma/s_j\\), a composition of both
hyperparameters per coordinate, which the map expresses exactly.

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
does not standardize and takes no such argument. Its columns are
grouping indicators and its penalty is a variance component with a
meaning of its own; weighting it by the size of the groups would change
the model.

## See also

[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
for the five penalties;
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
for the unpenalized block;
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
and [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md)
for the penalized structures;
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
and
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what a built one reports.

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20),
                 g = factor(rep(letters[1:4], 5)))

# The intercept is removed, so a four-level factor gives four columns.
built <- term_build(lasso(~ x1 + x2), dd)
term_coef_names(built)
#> [1] "lasso.x1" "lasso.x2"
term_coef_names(term_build(lasso(~ g), dd))
#> [1] "lasso.ga" "lasso.gb" "lasso.gc" "lasso.gd"

# Unless the intercept is all the formula has.
term_coef_names(term_build(ridge(~ 1), dd))
#> [1] "ridge.(Intercept)"

# The hyperparameters and the kink are the penalty's.
term_penalty(built)@params
#> [1] "lambda"
term_smooth(built)
#> [1] FALSE

# A matrix column of the data is the input that predicts. The name must
# resolve where the term is written as well as in the data.
R <- matrix(rnorm(60), 20, 3, dimnames = list(NULL, c("a", "b", "c")))
dd$R <- R
bm <- term_build(ridge(R), dd)
term_coef_names(bm)
#> [1] "ridge.a" "ridge.b" "ridge.c"
dim(term_predict(bm, dd[1:5, ]))
#> [1] 5 3

# Standardizing puts the spreads on the penalty's map, not on the design.
dd$x3 <- 1000 * dd$x2
bs <- term_build(lasso(~ x1 + x3, standardize = TRUE), dd)
term_penalty(bs)@map
#> 2 x 2 diagonal matrix of class "ddiMatrix"
#>           [,1]     [,2]
#> [1,] 0.7822996        .
#> [2,]         . 886.9947
apply(term_matrix(bs), 2, sd)
#>    lasso.x1    lasso.x3 
#>   0.7822996 886.9946868 

# An argument named after another penalty's hyperparameter is reported.
try(scad(~ x1, gamma = 1))
#> Error : 'scad' has no argument 'gamma'.
#>   Its hyperparameters are: lambda, a, each held by naming it and estimated when left NULL.
```
