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
  label = "random"
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
