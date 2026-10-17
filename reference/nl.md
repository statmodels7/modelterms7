# Nonlinear Parametric Term

Constructs a term for a contribution \\f(x; \theta)\\ that is nonlinear
in its own parameters and is given either as a formula in the covariates
and the parameters or as a function of both. Each parameter may carry a
link, and, when the term is given as a formula, may itself be modeled
with covariates.

## Usage

``` r
nl(
  fn,
  ...,
  params = NULL,
  x = NULL,
  links = NULL,
  subformulas = NULL,
  start = NULL,
  gradient = NULL,
  hessian = NULL,
  deriv3 = NULL,
  deriv4 = NULL,
  label = "nl"
)
```

## Arguments

- fn:

  A one-sided formula in the covariates and the parameters, or a
  function of `(x, theta)` vectorized in both.

- ...:

  Two-sided formulas whose left side names a parameter, one per
  parameter to be modeled with covariates, e.g. `theta1 ~ s(z)`. They
  apply to formula input only and have the same effect as naming the
  right-hand side in `subformulas`.

- params:

  The parameter names. They are required when `fn` is a function and
  rejected when `fn` is a formula, whose parameters are read from it.

- x:

  The covariate expression handed to a function `fn`, evaluated in the
  data. A formula does not use it.

- links:

  An optional named list of linkfunctions7 links, one per parameter.
  Parameters without one carry the identity.

- subformulas:

  An optional named list of one-sided formulas, the programmatic
  spelling of the formulas of `...`. It applies to formula input only.

- start:

  An optional named list of starting values for the parameters, on the
  parameter scale. A parameter that is not named starts at the inverse
  link of zero. A name that is not a parameter is rejected by
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- gradient, hessian, deriv3, deriv4:

  Optional functions `function(theta, data)` returning the derivatives
  of `fn` in its own parameters, each a named list of numeric vectors.
  Each may be supplied independently of the others; see the section on
  supplying the derivatives.

- label:

  A single non-empty string prefixed to the coefficient names.

## Value

An object of class
[`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md)
(a specification; see
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

While \\f\\ is differentiable the term is an ordinary additive one at
every point. The contribution is linearized as \$\$f(x;\theta(\beta))
\approx f(x;\theta(\beta_0)) + J(\beta_0)\\(\beta - \beta_0), \qquad J =
\frac{\partial f}{\partial \beta},\$\$ so the design block is the
Jacobian. The term differs from a linear one in that the block is
recomputed as the parameters move:
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
does that, and
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
reports the contribution itself, which is the second quantity that a
Gauss-Newton step needs.

### Two ways to give the function

A formula such as `~ theta1 * exp(theta2 * x)` is read symbolically: the
names in the formula that are not columns of the data are the
parameters, and the derivatives come from
[`stats::deriv()`](https://rdrr.io/r/stats/deriv.html) and
[`stats::D()`](https://rdrr.io/r/stats/deriv.html) where those succeed
and from a central difference where they do not. A function
`f(x, theta)`, vectorized in both, is treated as opaque: its derivatives
are the ones supplied in `gradient`, `hessian`, `deriv3` and `deriv4`
and are otherwise differenced, and its parameters must be named in
`params`.

The two routes differ in more than the derivatives. Modeling a parameter
with covariates means replacing \\\theta_j\\ by \\g_j^{-1}(Z\gamma_j)\\
inside \\f\\, which requires the place where \\\theta_j\\ enters. The
expression of a formula identifies this place and an opaque function
does not. `subformulas` is therefore available on the formula route
only, and is rejected on the other.

### Links and submodels

`links` carries each parameter to an unconstrained scale, so a rate
constrained positive is estimated as its logarithm and the optimizer
never proposes a negative one. A subformula develops a parameter as
\\\theta_j = g_j^{-1}(Z\gamma_j)\\ for a design \\Z\\ built from its
right-hand side, which gives a parameter that varies by group or with a
covariate; the coefficients are then the \\\gamma_j\\, and the Jacobian
carries the chain rule \\\partial f/\partial\theta_j \cdot (g_j^{-1})'
\cdot Z\\. Because the development acts on the unconstrained scale, the
parameter stays in its own set at every observation whatever the
coefficients are: a rate on the log link is positive for every row of
\\Z\\.

A subformula is written in `...` as a two-sided formula whose left side
names the parameter, `theta1 ~ z`, or programmatically as
`subformulas = list(theta1 = ~z)`; the two spellings are the same and a
parameter may carry only one. The right-hand side goes through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
so it takes any term of this package: `theta1 ~ ridge(~g)` is a
population value (the intercept) plus shrunken departures,
`theta1 ~ s(z)` lets the parameter move smoothly with a covariate, and
`theta1 ~ random(~1 | g)` is a random intercept on the parameter's
unconstrained scale. The penalties that the sub-terms carry are reported
through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
each entry named `parameter::subterm`, so a fitting layer estimates
their hyperparameters as it does those of any other term. A structural
term is rejected.

### A parameter that changes at a break-point

A subformula may carry a break-point term whose block is the Jacobian of
its contribution:
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
and
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
or [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
with `smoothed =`. Then \\g_j(\theta_j) = \eta_j(\gamma_j)\\ is no
longer linear in the coefficients: it is the sum of the sub-terms'
contributions. Under a log link,
`nl(~ a * exp(-r * x), r ~ jump(t, smoothed = sm), links = list(r = linkfunctions7::log_link()))`
with `sm` a
[`numericals7::smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.html)
is a rate whose logarithm steps at an estimated time. The block of \\f\\
carries the chain rule through the sub-term's own block, \\\partial
f/\partial\theta_j \cdot (g_j^{-1})' \cdot X_s(\gamma_s)\\, and its
derivatives in the coefficients add the sub-term's own
[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md),
[`term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.md)
and
[`term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.md)
where \\X_s\\ moves, so a fitting layer's exact outer derivatives reach
it unchanged. The position is confined to the sub-term's own interval,
as it is when the sub-term is used as a term by itself.

A sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
is rejected, naming `smoothed =`: its block is a working linearization
with a frozen weight and a position read off a product of two
coefficients, which is not a Jacobian and cannot sit inside another one.

The objective in a break-point's position has several optima, so where a
response is given,
[`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
tries each break-point at five interior quantiles of its covariate
(every combination of them for several break-points), takes damped
Gauss-Newton steps on the response from each and keeps the best fit. The
restarts that a fitting layer applies to a break-point term on its own
do not reach a break-point term nested here.

The fitted positions are read off the sub-term, which the term carries
at its fitted coefficients:
`seg_psi(term_components(built)$r$subs[[k]])`, with `k` the sub-term's
place in the subformula.

### Penalizing a parameter

A penalty is specified inside the subformula, by the sub-term that
carries it: `theta1 ~ lasso(~z1 + z2)` selects which covariates the
parameter depends on, and `theta1 ~ ridge(~g)` or
`theta1 ~ random(~1 | g)` shrinks per-group departures towards a
population value (the intercept, which stays unpenalized). Each sub-term
has its own hyperparameter, so that two parameters of a nonlinear
function can be penalized on scales of their own. `nl()` has no
`penalty` or `penalize` argument.

### Supplying the derivatives

The derivatives of \\f\\ in its own parameters determine the design
block (the first order), and a fitting layer reads the second and third
orders where the block moves with the coefficients. Each order is
obtained separately, by the first of three routes that applies: a
function supplied here, symbolic differentiation of the expression, and
one stencil applied to the highest analytic order below the one wanted,
or to the value of \\f\\ if every lower order is computed numerically. A
function given without derivatives, and a formula that
[`stats::deriv()`](https://rdrr.io/r/stats/deriv.html) and
[`stats::D()`](https://rdrr.io/r/stats/deriv.html) cannot differentiate,
therefore have every order differenced from the value.

`gradient`, `hessian`, `deriv3` and `deriv4` may be supplied
independently of one another. A supplied order is the starting point of
the stencil for the orders above it, so the third and fourth orders of a
function are one difference away from an exact second order instead of
from the value. A formula that
[`stats::D()`](https://rdrr.io/r/stats/deriv.html) can differentiate has
all four orders in closed form, whether or not a function is supplied
for some of them.

Each function takes `(theta, data)`: a named list of the parameters,
each of length `n` or 1, and the columns of the data as a list. It
returns a named list of numeric vectors of length 1 or `n`. The names
are the components of that order, keyed as distributions7 keys its own
derivative surfaces: the parameter names joined by `"_"`. Order two of
`c("a", "r")` is therefore `a_a`, `a_r`, `r_r`, and order three `a_a_a`,
`a_a_r`, `a_r_r`, `r_r_r`.

The names are normalized here, so `r_a` and `a_r` are the same component
and the order in which they are returned does not matter; the check
compares the set of names. A name that is not a component of this term,
a missing component, a repeated one or a component of a length other
than 1 or `n` is an error at
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md).
The term does not fall back to the numerical route in these cases.

The derivatives are taken in the parameters \\\theta\\ and not in the
coefficients. The term applies the chain rule onto the coefficients,
through each parameter's link and the design of its subformula.
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)
returns the derivatives of any order from one to four.

## See also

[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md),
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0, 3, length.out = 60))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(60, sd = 0.05)

# The term is an exponential decay, and a log link holds the rate positive.
spec <- nl(~ a * exp(-r * x),
           links = list(r = linkfunctions7::log_link()),
           start = list(a = 1, r = 1))
built <- term_build(spec, dd)
term_coef_names(built)
#> [1] "nl.a" "nl.r"
dim(term_matrix(built))
#> [1] 60  2

# The amplitude is developed by group: a population value plus shrunken
# departures, whose hyperparameter a fitting layer estimates.
dd$g <- factor(rep(c("u", "v"), 30))
sub <- term_build(nl(~ a * exp(-r * x), a ~ ridge(~g),
                     start = list(a = 1, r = 1)), dd)
vapply(term_penalties(sub), function(e) e$name, character(1))
#> [1] "a::ridge(~g)"

# The first two orders are written out by hand.
g <- function(theta, data) {
  e <- exp(-theta$r * data$x)
  list(a = e, r = -theta$a * data$x * e)
}
h <- function(theta, data) {
  e <- exp(-theta$r * data$x)
  list(a_a = 0 * e, r_r = theta$a * data$x^2 * e, a_r = -data$x * e)
}
ex <- term_build(nl(~ a * exp(-r * x), gradient = g, hessian = h,
                    start = list(a = 1, r = 1)), dd)
names(nl_fderiv(ex, order = 1))
#> [1] "a" "r"

# for a function, the orders that are not supplied are differenced from
# the highest one supplied
f <- function(x, theta) theta$a * exp(-theta$r * x)
exf <- term_build(nl(f, params = c("a", "r"), x = x, gradient = g,
                     hessian = h, start = list(a = 1, r = 1)), dd)
max(abs(nl_fderiv(exf, order = 3)$r_r_r + dd$x^3 * exp(-dd$x)))
#> [1] 9.263701e-13

# The model is fitted to data simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(8)
  fd <- data.frame(x = sort(runif(150, 0, 5)))
  fd$y <- 3 * exp(-0.8 * fd$x) + rnorm(150, sd = 0.1)
  ft <- statmodels7::statmod(
    y ~ 0 + nl(~ a * exp(-r * x), start = list(a = 1, r = 1)),
    distributions7::gaussian1_distrib(), fd)
  # The true values are a = 3 and r = 0.8.
  round(coef(ft)$mu, 2)
}
#> nl.a nl.r 
#> 2.98 0.80 
```
