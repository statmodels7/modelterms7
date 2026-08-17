# Penalized Smooth of One Covariate

A smooth function of a covariate, expanded in a basis7 basis and
penalized for roughness. The default construction is a cubic B-spline
basis under the Demmler-Reinsch reparametrization, which separates the
linear effect from the nonlinear deviation and turns the roughness
penalty into the identity on the deviation.

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
  sparse = NULL
)
```

## Arguments

- x:

  The covariate, an expression evaluated in the data.

- by:

  An optional factor or numeric variable; see Details.

- k:

  The basis dimension before reparametrization. Defaults to 10.

- degree:

  The spline degree. Defaults to 3, a cubic spline.

- basis:

  An optional basis7 basis to use in place of the default B-spline; its
  range is taken as given.

- linear:

  Whether the linear effect is carried in the block, and left
  unpenalized there. Defaults to `TRUE`.

- label:

  A single non-empty string prefixed to the coefficient names. Defaults
  to a name built from the covariate.

- lambda:

  The smoothing parameter, held at the value given and ESTIMATED when
  left `NULL`, which is the default. An anisotropic tensor product
  carries one per margin, so a vector of that length, or a named one
  holding some of them.

- sparse:

  Whether the block is built as a `dgCMatrix`. `NULL`, the default,
  settles it at build from the size of the block. Only a FACTOR `by`
  admits it, each row sitting in the block of its own level; without one
  an explicit `TRUE` is refused rather than ignored. See Details.

## Value

An object of class
[`SmoothTerm`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The block has one column for the linear effect, centered and scaled,
followed by the reparametrized basis. That ordering is what the penalty
reads: it is the quadratic penalty of \\\mathrm{diag}(0, 1, \dots, 1)\\,
rank deficient by exactly one, so the linear effect is unpenalized and
the deviation is shrunk towards zero. As the smoothing parameter grows
the fit approaches a straight line rather than a constant, and
[`edf`](https://statmodels7.github.io/modelterms7/reference/edf.md)
falls towards one.

The Demmler-Reinsch construction (demmler1975, as used for effect
selection by bach2024) is empirical: it takes the inner product at the
observed covariate values, so it is built when the term is, and the
transform is stored in the blueprint. Prediction at new points is the
parent basis evaluated there and multiplied by the same transform, so
the separation of the linear from the nonlinear part holds on new data
as it does on old.

### Varying the smooth by another variable

With `by` a factor the term carries one smooth per level, the block
being the smooth multiplied by each level's indicator and the penalty
the same matrix repeated blockwise, so one smoothing parameter governs
every level. With `by` numeric the term is a varying-coefficient one:
the smooth multiplies that variable, and the fitted function is the
coefficient of `by` as it changes with the covariate.

A factor `by` is where a smooth's block can be SPARSE: each row sits in
the block of its own level and nowhere else, a density of \\1/m\\.
`sparse = TRUE` builds it that way rather than building the dense matrix
and compressing it – measured at 2000 rows, \\k = 10\\ and 200 levels,
0.35 MB against 28.93 MB, the numbers identical. `sparse = NULL`, the
default, settles it at build from the size of the block, the dense form
holding \\n m k\\ cells against \\n k\\ non-zeros; see
[`.resolve_sparse`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md)
for the threshold and what it was measured against. An explicit `TRUE`
is refused without a factor `by`, where there would be nothing to build
on: the basis is dense by construction, the Demmler-Reinsch rotation
making it so, and a numeric `by` merely multiplies it.

The block alone is sparse. The PENALTY of a factor `by` is the same
matrix repeated blockwise, and penalties7 returns it dense – 25.92 MB at
those sizes, at a density of 0.0005 – which is a property of that
package's contract rather than of this construction.

## References

Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
smoothing. *Numerische Mathematik*, 24, 375–382.

Bach, P. and Klein, N. (2024). Bayesian effect selection in additive
models with an application to time-to-event data.

## See also

[`te`](https://statmodels7.github.io/modelterms7/reference/te.md),
[`random`](https://statmodels7.github.io/modelterms7/reference/random.md),
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(80)))
dd$y <- sin(2 * pi * dd$x) + rnorm(80, sd = 0.2)
built <- term_build(s(x, k = 8), dd)
term_npar(built)
#> [1] 7
term_penalty(built)@params
#> [1] "lambda"
```
