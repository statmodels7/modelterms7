# Penalized Smooth of One Covariate

A smooth function of a covariate, expanded in a basis7 basis and
penalized for roughness. The default construction is a cubic B-spline
basis under the Demmler-Reinsch reparametrization, which separates the
linear effect from the nonlinear deviation and turns the roughness
penalty into the identity on the deviation.

## Usage

``` r
s(x, by = NULL, k = 10, degree = 3, basis = NULL, linear = TRUE, label = NULL)
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
