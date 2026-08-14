# Grouped Random-Effect Term

Random intercepts and slopes for a grouping factor: `random(~ 1 | g)`
builds one coefficient per level of `g`, and `random(~ x | g)` one
intercept and one slope per level, with the distribution of the effects
attached as the penalty on those coefficients – which is what a random
effect is under penalized likelihood.

## Usage

``` r
random(
  formula,
  correlated = TRUE,
  precision = NULL,
  distrib = NULL,
  kinks = numeric(0),
  label = "random",
  hyper = NULL
)
```

## Arguments

- formula:

  A bar formula, `~ lhs | g`, with `g` evaluating to the grouping
  variable in the data.

- correlated:

  Logical; whether the default Gaussian lets the within-group effects
  correlate. Ignored when `precision` or `distrib` is given.

- precision:

  A parameters7 matrix parameter of the within-group dimension, or
  `NULL`.

- distrib:

  A univariate distributions7 object for the effects, or `NULL`;
  exclusive with `precision`.

- kinks:

  The kink set of `distrib`'s log-density in its argument, passed to
  [`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.html).

- label:

  A single non-empty string prefixed to the coefficient names.

- hyper:

  The hyperparameters of the effects' distribution to HOLD, as a named
  vector or list; those not named are estimated. Which names there are
  depends on what the term was given: `c(sigma = 0.4)` for the default,
  the free names of the structure for a `precision`, and the
  distribution's own parameters for a `distrib`. A name the penalty does
  not carry is reported when the term is built, which is where the
  penalty first exists.

## Value

An object of class
[`RandomTerm`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The left side of the bar is an ordinary one-sided formula for the
within-group design, with the usual intercept convention: `~ x | g`
carries an intercept and a slope per group and `~ 0 + x | g` the slope
alone. The block interacts that design with the group indicators,
ordered group by group, so the coefficients of one group are adjacent.

Three choices of effect distribution are available. By default the
effects are Gaussian with one covariance shared by every group:
unstructured over the within-group coefficients when `correlated = TRUE`
(a
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
chart, so intercepts and slopes may correlate), diagonal when
`correlated = FALSE` (independent effects, one variance per within-group
column), and a single scale when the bar's left side is the intercept
alone. With `precision`, a parameters7 matrix parameter of the
within-group dimension replaces that default (an AR(1) over ordered
within-group columns, a compound symmetry); it enters as the per-group
precision, replicated across groups by
[`kron_identity`](https://statmodels7.github.io/parameters7/reference/kron_identity.html)
into
[`structured_penalty`](https://statmodels7.github.io/penalties7/reference/structured_penalty.html),
and its hyperparameters are the structure's free values. With `distrib`,
a univariate distributions7 object – holding its own location, typically
through
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.html)
at zero – is applied coordinatewise to every effect through
[`distrib_penalty`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.html).
A distribution used as a penalty gives joint-mode (penalized likelihood)
estimation of the effects; the marginal likelihood, which integrates
them out, is not provided here.

Prediction maps new data onto the levels seen at build time; a level the
term has not seen is rejected.

A random effect is not standardized, and there is no `standardize`
argument to ask for it with; passing one is an error. The term is a
ridge where the effects are Gaussian and independent (see below), but
its columns are grouping indicators rather than measured covariates, and
its hyperparameter is a variance component with a meaning of its own.
Dividing each coefficient by the spread of its indicator would weight
the effects by the sizes of the groups, which changes the model rather
than the scale its hyperparameter is read on.

## The block and its penalty

With \\m\\ levels and a within-group design \\Z_i\\ of \\d\\ columns,
the block is the interaction of that design with the group indicators,
ordered group by group,

\$\$Z = \operatorname{diag}(Z_1, \dots, Z_m), \qquad b = (b_1', \dots,
b_m')',\$\$

so the \\d\\ coefficients of one group occupy adjacent positions. Under
the default Gaussian the penalty is the negative log-density of \\b\\ at
a precision replicated across groups,

\$\$\rho(b; \eta) = \tfrac{1}{2} b' \left(I_m \otimes
\Omega_d(\eta)\right) b - \tfrac{m}{2}\log\lvert \Omega_d(\eta)\rvert +
\tfrac{md}{2}\log(2\pi),\$\$

whose hyperparameters are the free values \\\eta\\ of the within-group
parameter. Minimizing the penalized least squares in \\(\beta, b)\\ is
the mixed-model equation at the variance ratio the precision encodes, so
the minimizer is the best linear unbiased predictor; \\\Omega_d =
I_d/\sigma_b^2\\ recovers the ordinary random intercept, where the
penalty is exactly the ridge.

## References

Laird, N. M. and Ware, J. H. (1982). Random-effects models for
longitudinal data. *Biometrics* 38, 963-974.

## See also

[`s`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te`](https://statmodels7.github.io/modelterms7/reference/te.md),
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md)

## Examples

``` r
dd <- data.frame(y = rnorm(9), x = rnorm(9),
                 g = factor(rep(c("a", "b", "c"), 3)))
built <- term_build(random(~ x | g), dd)
term_coef_names(built)
#> [1] "random.a.(Intercept)" "random.a.x"           "random.b.(Intercept)"
#> [4] "random.b.x"           "random.c.(Intercept)" "random.c.x"          
term_penalty(built)@params
#> [1] "log_L1" "log_L2" "L2.1"  
```
