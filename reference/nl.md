# Nonlinear Parametric Term

A contribution \\f(x; \theta)\\ that is nonlinear in its own parameters,
given either as a formula in the covariates and the parameters or as a
function of both. Each parameter may carry a link, and, when the term is
given as a formula, may itself be modeled with covariates.

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
  label = "nl"
)
```

## Arguments

- fn:

  A one-sided formula in the covariates and the parameters, or a
  function of `(x, theta)` vectorized in both.

- ...:

  Two-sided formulas whose left side names a parameter, one per
  parameter to be modeled with covariates, e.g. `theta1 ~ s(z)`. Formula
  input only; the same as naming the right-hand side in `subformulas`.

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

  An optional named list of one-sided formulas, the programmatic
  spelling of the formulas of `...`. Formula input only.

- start:

  An optional named list of starting values for the parameters, on the
  parameter scale. Defaults to the inverse link at zero.

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
[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
so it takes any term of this package: `theta1 ~ ridge(g)` is a
population value (the intercept) plus shrunken departures,
`theta1 ~ s(z)` lets the parameter move smoothly with a covariate, and
`theta1 ~ random(~1 | g)` is a random intercept on the parameter's
unconstrained scale. The penalties the sub-terms carry are reported
through
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
under the key `parameter::subterm`, so a fitting layer estimates their
hyperparameters as it does any other term's. A structural term, and a
term whose block moves with its coefficients, are rejected: a
parameter's submodel must be a fixed design.

### Penalizing a parameter

A penalty is asked for inside the subformula, where the sub-term that
carries it declares it: `theta1 ~ lasso(~z1 + z2)` selects which
covariates the parameter depends on, and `theta1 ~ ridge(~g)` or
`theta1 ~ random(~1 | g)` shrinks per-group departures towards a
population value (the intercept, which stays unpenalized). Each sub-term
brings its own hyperparameter, which is what two parameters of a
nonlinear function want, being on scales of their own. Earlier releases
carried `penalty` and `penalize` arguments over the parameters'
coefficients; both are gone, the sub-terms covering the cases that
matter with the hyperparameters in the right place.

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

# the amplitude developed by group: a population value plus shrunken
# departures, whose hyperparameter a fitting layer estimates
dd$g <- factor(rep(c("u", "v"), 30))
sub <- term_build(nl(~ a * exp(-r * x), a ~ ridge(~g),
                     start = list(a = 1, r = 1)), dd)
vapply(term_penalties(sub), function(e) e$name, character(1))
#> [1] "a::ridge(~g)"
```
