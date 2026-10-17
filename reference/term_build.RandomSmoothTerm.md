# Build a Smooth Random Effect

Builds the smooth on the left side of the bar, separates its null space
from its penalized coordinates, interacts the resulting within-group
design with the group indicators, and attaches one prior to the null
space and one to the penalized coordinates.

## Arguments

- term:

  An unbuilt or built
  [`RandomSmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomSmoothTerm.md).

- data:

  A data frame carrying the grouping variable and the smooth's
  covariates.

- ...:

  Unused.

## Value

The term with `X` (a `dgCMatrix`), `coef_names`, `blueprint` and
`penalty` filled.

## The within-group design

The design that the effect of one group multiplies is the constant
(unless the left side carries `0`), then the null space of the smooth's
penalty without the constant, which every smoother removes, then the
penalized coordinates, as
[`.random_smooth_split()`](https://statmodels7.github.io/modelterms7/reference/dot-random_smooth_split.md)
defines them. The null space is whatever the smoother's operator leaves
unpenalized: a straight line for a second derivative, a sine and a
cosine for the harmonic operator, an exponential for \\D^2 + aD\\,
nothing for a cyclic smoother.

The argument `null_space` of the smoother determines what the null space
contributes. Under `"keep"` its columns are there and carry the
null-space prior together with the constant; under `"shrink"` they are
penalized with the rest and belong to the penalized coordinates; under
`"drop"` they are absent.

## The two priors

The constant and the null-space columns, \\d_0\\ of them, have the prior
that
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
gives \\d_0\\ within-group columns: a centered multivariate Gaussian on
an unstructured covariance, or one standard deviation per column under
`correlated = FALSE`, or `distrib` where it is given. A full covariance
makes the model invariant to the origin and the scale of the covariate,
and to the phase of a periodic one.

The penalized coordinates of every group have the prior of the smooth's
own penalty. With one component that is an independent centered Gaussian
with one standard deviation, because the coordinates have been whitened,
and its hyperparameter is `sigma.smooth`. With several components it is
the sum \\\sum_k \lambda_k (I_m \otimes S_k)\\, with hyperparameters
`lambda1.smooth`, `lambda2.smooth` and so on. The groups share these
hyperparameters, so all the curves have the same roughness.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md),
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40), g = factor(rep(letters[1:4], 10)))
b <- term_build(random(~ s(x, smoother = basis7::bspline_smooth(k = 5)) | g),
                dd)
term_group(b)$names
#> [1] "(Intercept)" "lin"         "z1"          "z2"          "z3"         
vapply(term_penalties(b), function(e) paste(e$penalty@params, collapse = ","),
       character(1))
#> [1] "sigma_log_L1,sigma_log_L2,sigma_L2.1"
#> [2] "sigma"                               
```
