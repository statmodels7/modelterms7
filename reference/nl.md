# Nonlinear Parametric Term

A contribution \\f(x; \theta)\\ that is nonlinear in its own parameters,
given either as a formula in the covariates and the parameters or as a
function of both. Each parameter may carry a link, and, when the term is
given as a formula, may itself be modeled with covariates.

## Usage

``` r
nl(
  fn,
  params = NULL,
  x = NULL,
  links = NULL,
  subformulas = NULL,
  start = NULL,
  penalty = c("none", "lasso", "ridge"),
  penalize = NULL,
  label = "nl"
)
```

## Arguments

- fn:

  A one-sided formula in the covariates and the parameters, or a
  function of `(x, theta)` vectorized in both.

- params:

  The parameter names. Required when `fn` is a function; inferred from
  the formula otherwise.

- x:

  The covariate expression handed to a function `fn`, evaluated in the
  data. Unused for a formula.

- links:

  An optional named list of linkfunctions7 links, one per parameter.
  Parameters without one carry the identity.

- subformulas:

  An optional named list of one-sided formulas, one per parameter to be
  modeled with covariates. Formula input only.

- start:

  An optional named list of starting values for the parameters, on the
  parameter scale. Defaults to the inverse link at zero.

- penalty:

  One of `"none"` (default), `"lasso"` or `"ridge"`, applied to the
  coefficients of the parameters `penalize` names.

- penalize:

  The parameters the penalty reaches, as a character vector. Defaults to
  all of them.

- label:

  A single non-empty string prefixed to the coefficient names.

## Value

An object of class
[`NlTerm`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

While \\f\\ is differentiable the term is an ordinary additive one at
every point: the contribution is linearized as \$\$f(x;\theta(\beta))
\approx f(x;\theta(\beta_0)) + J(\beta_0)\\(\beta - \beta_0), \qquad J =
\frac{\partial f}{\partial \beta},\$\$ so the design block is the
Jacobian, and the only thing that distinguishes the term from a linear
one is that the block is refreshed as the parameters move.
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
does that, and
[`term_value`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
reports the contribution itself, which a Gauss-Newton step needs beside
the Jacobian.

### Two ways to give the function, with different reach

A **formula** such as `~ theta1 * exp(theta2 * x)` is read symbolically:
the names it uses that are not columns of the data are the parameters,
and the derivatives come from
[`deriv`](https://rdrr.io/r/stats/deriv.html) where that succeeds and
from a central difference where it does not. A **function**
`f(x, theta)`, vectorized in both, is treated as opaque: its derivatives
are always differenced, and its parameters must be named in `params`.

The difference is not only in the derivatives. Modeling a parameter with
covariates means replacing \\\theta_j\\ by \\g_j^{-1}(Z\gamma_j)\\
inside \\f\\, which requires knowing where \\\theta_j\\ enters; a
formula says so and an opaque function does not. `subformulas` is
therefore available on the formula route only, and is rejected on the
other.

### Links and submodels

`links` carries each parameter to an unconstrained scale, so a rate
constrained positive is estimated as its logarithm and the optimizer
never proposes a negative one. `subformulas` develops a parameter as
\\\theta_j = g_j^{-1}(Z\gamma_j)\\ for a design \\Z\\ built from a
one-sided formula, which gives a parameter that varies by group or with
a covariate; the coefficients are then the \\\gamma_j\\, and the
Jacobian carries the chain rule \\\partial f/\partial\theta_j \cdot
(g_j^{-1})' \cdot Z\\.

### Penalizing a parameter

`penalty` attaches a penalty to the coefficients of the parameters
`penalize` names, one penalty per parameter and so one hyperparameter
each, since two parameters of a nonlinear function are on scales of
their own and have no reason to share one. They are declared through
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
which names the coefficients each covers; the parameters left out are
unpenalized.

What is shrunk is the coefficient, not the parameter, so with a link the
target is \\g_j^{-1}(0)\\ rather than zero: a rate carried by a log link
is shrunk towards one. Where the parameter carries a subformula, the
whole vector \\\gamma_j\\ is covered, so a lasso there selects which
covariates a parameter depends on. A subformula of the form `~ g` with a
factor `g` is the population-and-deviations pattern: the intercept is
the population value and the remaining columns are deviations, penalized
as the coordinates they are.

## See also

[`s`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te`](https://statmodels7.github.io/modelterms7/reference/te.md),
[`random`](https://statmodels7.github.io/modelterms7/reference/random.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0, 3, length.out = 60))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(60, sd = 0.05)

# an exponential decay, with the rate held positive by a log link
spec <- nl(~ a * exp(-r * x),
           links = list(r = linkfunctions7::log_link()),
           start = list(a = 1, r = 1))
built <- term_build(spec, dd)
term_coef_names(built)
#> [1] "nl.a" "nl.r"
dim(term_matrix(built))
#> [1] 60  2
```
