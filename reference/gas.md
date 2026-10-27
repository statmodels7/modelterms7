# Score-Driven Dynamics

`gas()` constructs a generalized autoregressive score component (Creal
et al., 2013; Harvey, 2013): a level \\f_t\\ that is added to the linear
predictor and driven by the score of the observation density at the
previous times, \$\$f_t = \omega + \sum\_{i=1}^{p} \kappa_i s\_{t-i} +
\sum\_{j=1}^{q} \xi_j f\_{t-j},\$\$ with \\s_t = \partial \ell_t /
\partial \eta_t\\ the derivative of the log-likelihood contribution at
time \\t\\ with respect to the predictor at which it is evaluated,
multiplied by a power of the expected information when `scaling` is not
zero.

## Usage

``` r
gas(
  p = 1,
  q = 1,
  ...,
  by = NULL,
  time = NULL,
  scaling = 0,
  p_cross = NULL,
  links = NULL,
  label = "gas"
)
```

## Arguments

- p:

  The number of score lags, a finite integer of at least 1. Defaults to
  1.

- q:

  The number of autoregressive lags, a finite integer of at least 0.
  Defaults to 1.

- ...:

  Two-sided formulas whose left side names a parameter, one per
  parameter to be developed with covariates, e.g. `kappa1 ~ s(x)`; see
  the section above. The left side may also name a cross loading,
  `kappa<i>_<parameter>`, whose lag `i` the term declares in `p_cross`,
  e.g. `kappa1_mu ~ 1 + random(~ 1 | g)`; the cross loading is developed
  on the identity link unless `links` names another. A cross loading
  that varies by group is identified only with long series or many
  groups (with a few hundred observations per group). A subformula with
  no column, such as `kappa2_mu ~ 0`, is rejected: a cross loading is
  restricted to zero by leaving its lag out of `p_cross`. A formula
  given in the place of `p` or `q` is rejected.

- by:

  An optional grouping variable, evaluated in the data; each group is
  filtered independently, from its own starting level. A formula here is
  the shorthand that gives the same subformula to every parameter of the
  filter (not to the cross loadings), and the term then does not imply
  any grouping.

- time:

  An optional ordering variable, evaluated in the data.

- scaling:

  A single finite number \\d\\: the score is multiplied by the expected
  information to the power \\-d\\. Defaults to 0, the score as it is;
  see the section on the scaling of the score.

- p_cross:

  An optional named vector of integers of at least 0, one per other
  filtered parameter whose score the level of this term loads, giving
  the number of lags of that score, e.g. `p_cross = c(sigma = 2)` for
  the loadings `kappa1_sigma` and `kappa2_sigma`. Each named parameter
  must carry a `gas()` term of its own; see the section on several
  filtered parameters. `NULL`, the default, gives a term that loads its
  own score only, and so does a count of 0.

- links:

  An optional named list of linkfunctions7 links over the parameters of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
  overriding the defaults described above. Only the parameters of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  before any development, and the cross loadings `kappa<i>_<parameter>`,
  can be named. The coefficients of a development are unconstrained, and
  the link of the developed parameter is applied inside the development.

- label:

  A single non-empty string naming the term.

## Value

An object of class
[`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md)
(a specification; see
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The term adds no columns. The predictor at one time depends on the data
at the previous ones, so the contribution cannot be written as a block,
and
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
runs the recursion instead. For this reason the term belongs to the
structural branch.

The score of the distribution that the model carries drives the
recursion, so the same term is a GARCH-like volatility model when it
enters the scale of a Gaussian, a dynamic count model when it enters the
mean of a Poisson, and a robust location filter when it enters a Student
t. A heavy-tailed score is bounded in the observation, so an outlier
moves the level by a bounded amount instead of in proportion to its
size.

### The parameters and their chart

The parameters are the level \\\omega\\, the score loadings \\\kappa_1,
\dots, \kappa_p\\, and the persistence. Each is estimated on the
unconstrained scale of a link, and `links` overrides any of them; the
defaults are the following.

The level carries the identity link because it is unconstrained. The
loadings carry the log link: a positive loading responds in the
direction of the score, which is the case that the score-driven
literature writes, and positivity is then structural. A subformula
develops the loading on the log scale, so every group and every
observation has a positive loading. A loading that must be free in sign
takes the identity link, which is set with
`links = list(kappa1 = linkfunctions7::identity_link())`.

The persistence is carried by partial autocorrelations and not by the
coefficients \\\xi_j\\ directly. The stationary region of an
autoregression is not a box, so a collection of scalar links cannot
represent it, whereas the partial autocorrelations each range over
\\(-1, 1)\\ independently and the Levinson-Durbin recursion maps them
bijectively onto the coefficients. At \\q = 1\\ the two coincide. Each
coordinate is named after the chart on which it is estimated, `pacf1`
and so on, following the convention of parameters7.

With the default links, a parameter developed over groups or
observations stays in its own set. The development acts on the
unconstrained scale, so a loading on the log link is positive and a
persistence on the rhobit link is stationary at every observation.

### Several filtered parameters

A `gas()` term in the formula of each of several parameters, as in
`y ~ gas(p = 1, q = 1, time = t) | sigma ~ gas(p = 1, q = 1, time = t)`,
gives one recursion for the vector of their levels. The score of each
filtered parameter is read at the levels of all the others at the same
time, so the levels cannot be filtered one after the other.
[`gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.md)
joins the terms, and a fitting layer reads them as one structural term.
The terms must walk the same series (the same `by` and `time`).

The scaling is joint. With \\G\\ the filtered parameters, the driving
vector is \$\$u_t = \mathcal{I}\_{GG,t}^{-d}\\ \frac{\partial
\ell_t}{\partial \eta\_{G,t}},\$\$ where \\\mathcal{I}\_{GG,t}\\ is the
expected information of the filtered parameters with the other
parameters held. Every term carries the same `scaling`, and with more
than one term it is 0, 1/2 or 1. For 1/2 the power is the symmetric
square root, which does not depend on the order of the parameters. Where
the information is diagonal, as for the mean and the scale of a
gaussian, the joint scaling equals the scaling of each parameter by its
own information.

By default each level loads its own score only, so the matrices of
loadings are diagonal. The argument `p_cross` of a term names the other
filtered parameters whose scores its level loads as well, with the
number of lags of each: in
`y ~ gas(p = 1, q = 1, p_cross = c(sigma = 1)) | sigma ~ gas(p = 1, q = 1)`
the level of `mu` responds to the score of `sigma` at lag one, and the
level of `sigma` responds to its own score only, so the matrix of
loadings at lag one is upper triangular. A lag that `p_cross` leaves out
is a loading restricted to zero, and the number of lags of a partner may
exceed the term's own `p`. The parameter named in `p_cross` must carry a
`gas()` term of its own. A cross loading rides the identity link,
because its sign is free. For a gaussian, the loading of the score of
the mean in the level of the log-scale is the asymmetric response of the
volatility to the sign of a surprise, which is the term in \\z\_{t-1}\\
of an EGARCH model. Each level keeps its own persistence, so the
autoregressive part is stationary when the part of each level is.

### The scaling of the score

The recursion is driven by the score multiplied by a power of the
expected information of the same predictor, \$\$s_t =
\mathcal{I}\_t^{-d}\\\frac{\partial \ell_t}{\partial \eta_t}, \qquad
\mathcal{I}\_t = -\mathrm{E}\\\left\[ \frac{\partial^2 \ell_t}{\partial
\eta_t^2}\right\],\$\$ where \\d\\ is the argument `scaling` and
\\\mathcal{I}\_t\\ is evaluated at the predictor and at the other
parameters of observation \\t\\. The level is scalar, so
\\\mathcal{I}\_t\\ is a number and every real \\d\\ defines a valid
model. The literature uses three values (Creal et al., 2013). With \\d =
0\\, the default, the score is used as it is. With \\d = 1/2\\ the
scaled score has unit conditional variance. With \\d = 1\\ the update is
the score premultiplied by the inverse information, a step of Fisher
scoring, and to first order it does not depend on the link of the
parameter. The value \\d = 1\\ with the identity link gives three known
models exactly: GARCH(1, 1) for the variance of a gaussian, the
INGARCH(1, 1) model for the mean of a Poisson and the ACD(1, 1) model
for the mean of an exponential.

Wherever \\\mathcal{I}\_t\\ varies with the observation, the loadings
\\\kappa_i\\ have different units for different values of \\d\\. Where
it is constant, as for the log-scale of a gaussian, changing \\d\\ only
rescales the loadings and leaves the likelihood unchanged. A scaling
other than zero requires a univariate family whose expected information
has an analytic second derivative
([`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html)),
because the derivatives of the filter reach the fourth derivative of
\\\mathcal{I}\_t\\, which is read from the derivatives of the expected
information up to
[`distributions7::distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.html).

\\\mathcal{I}\_t\\ is part of the model and is read in closed form at
the predictor of every step of the recursion. Every family of
distributions7 carries the closed form. A family written with its
density alone supplies it through methods for
[`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html),
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
and
[`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html);
without them only `scaling = 0` is available. The two substitutes of
distributions7 for a missing closed form do not apply here: the outer
product of the scores
([`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html))
reads \\\mathcal{I}\_t\\ at the data as \\(\partial \ell_t/\partial
\eta_t)^2\\, which gives \\s_t = (\partial \ell_t/\partial \eta_t)\\
\|\partial \ell_t/\partial \eta_t\|^{-2d}\\ and so defines another
model, and a quadrature
([`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html))
would be evaluated once per step and per derivative order, in sequence.

### The level beside an intercept

The level \\\omega\\ adds a constant to the recursion, so an intercept
\\c\\ in the same equation spans the same direction. The two are exactly
confounded: shifting \\c\\ by \\k\\ and \\\omega\\ by \\-k(1 - \sum_j
\xi_j)\\ leaves every predictor unchanged. A fitting layer therefore
estimates at most one of them. In statmodels7 the intercept is kept and
\\\omega\\ is held at zero, which a summary reports as `omega (held)`;
written without the intercept, as `0 + gas(...)`, the level is estimated
instead. The same holds for any set of columns whose span contains the
constant, such as a factor coded without an intercept.

The two spellings are one model with the same log-likelihood, loadings
and persistence, and their constants measure different quantities. The
score has zero conditional mean, so the level fluctuates around the
fixed point of the recursion, \$\$\mathrm{E}\[f_t\] = \frac{\omega}{1 -
\sum\_{j=1}^{q} \xi_j},\$\$ and the intercept estimates this quantity
when \\\omega\\ is held. Therefore \$\$\eta_t = c + f_t \\\text{with}\\
\omega = 0 \qquad\text{and}\qquad \eta_t = f_t \\\text{with}\\ c =
\frac{\omega}{1 - \sum_j \xi_j}\$\$ describe the same predictor. With a
persistence close to one the factor \\1/(1 - \sum_j \xi_j)\\ is large,
so \\\omega\\ is small, and its standard error is not comparable with
the intercept's: a change in \\\omega\\ moves the level by that factor.
The identity is for a scalar level; a level developed by a subformula is
confounded with the part of its span that the equation's design shares.

### Groups and time

`by` filters each group independently, as a panel of short series needs,
and `time` gives the order within a group. Only the order of the times
matters: a gap between two times counts as one step of the recursion.
Without `time` the rows are taken in the order in which they appear.

### A parameter developed with covariates

A two-sided formula in `...` whose left side names a parameter develops
it as \\\psi\_{j,t} = g_j^{-1}(z_t^\top\gamma_j)\\, where the design
\\Z\\ is built from the right-hand side through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
so the design takes the additive terms of the package whose block does
not move with their coefficients:

    gas(p = 1, q = 1, omega ~ ridge(~g), kappa1 ~ s(x),
        pacf1 ~ random(~1 | id), by = id)

The development acts on the unconstrained scale of the parameter's own
link, and therefore every per-observation value stays inside the
parameter's own set whatever the coefficients are: a loading on the log
link is positive at every observation, a persistence on the rhobit chart
is inside \\(-1, 1)\\ at every observation, and at \\q = 1\\ that bounds
the recursion's growth step by step. The coefficients \\\gamma_j\\ are
the term's parameters, unconstrained and on the identity link; the
penalties that the sub-terms carry are reported through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
under the key `parameter::subterm`. A parameter that varies by
observation changes the recursion itself, \$\$f_t = \omega_t + \sum_i
a\_{i,t}\\ s\_{t-i} + \sum_j b\_{j,t}\\ f\_{t-j},\$\$ with \\b_t\\ from
the Levinson-Durbin map of that observation's partial autocorrelations,
and the filter, its derivative, the reverse pass and the curvature all
run the general recursion.

`by = ~f` (a formula, where a grouping variable is a bare symbol) is the
shorthand giving the same subformula to every parameter; mixing it with
per-parameter formulas is an error. A structural term, and a term whose
block moves with its own coefficients
([`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
and
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)),
are rejected inside a subformula.

### A population value and a departure per group

The panel case is one subformula:
`gas(omega ~ random(~1 | id), by = id)` is a population value (the
intercept of the development) plus one unconstrained departure per
group, shrunk by the random intercept's own ridge, whose hyperparameter
a fitting layer estimates. `lasso(...)` in the subformula sets the
departures of the groups that do not need one exactly to zero.

The departures act on the unconstrained scale of the parameter's chart,
so a persistence stays inside \\(-1, 1)\\ whatever the departure is. The
penalty on the departures also identifies them: a population value and
\\m\\ departures are \\m+1\\ numbers describing \\m\\ group values, so
without the penalty the likelihood is flat along one direction per
developed parameter. This is the parametrization of a random effect and
it is identified the same way, there by a variance component, here by
the penalty. An unpenalized development over group indicators is
therefore for reading a filter at given parameters rather than for
fitting one.

The arguments `deviations` and `penalty` are not accepted and are
rejected by name; the subformula takes their place.

## References

Creal, D., Koopman, S. J. and Lucas, A. (2013). Generalized
autoregressive score models with applications. *Journal of Applied
Econometrics*, 28(5), 777–795.

Harvey, A. C. (2013). *Dynamic Models for Volatility and Heavy Tails*.
Cambridge University Press.

## See also

[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
for the other dynamic term,
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
for the recursion,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
and
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the parameters and their charts,
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
for what a fitted one reports.

## Examples

``` r
# The term has a level, one loading and one persistence coordinate.
term_params(gas(p = 1, q = 1))
#> [1] "omega"  "kappa1" "pacf1" 
term_params(gas(p = 2, q = 2))
#> [1] "omega"  "kappa1" "kappa2" "pacf1"  "pacf2" 

# The loading rides a log chart so that it stays positive, and the
# persistence a rhobit so that the filter stays stationary.
vapply(term_links(gas(p = 1, q = 2)), function(l) l@link_name,
       character(1))
#>      omega     kappa1      pacf1      pacf2 
#> "identity"      "log"   "rhobit"   "rhobit" 

# Running the filter: the term supplies the recursion and the caller
# supplies the density's score and curvature. Here the predictor is a
# Gaussian mean.
set.seed(1)
dd <- data.frame(t = 1:60, y = c(rnorm(30), rnorm(30, 3)))
b <- term_build(gas(p = 1, q = 1, time = t), dd)
out <- term_filter(b, eta = rep(0, 60), y = dd$y,
                   score = function(e, i) dd$y[i] - e,
                   curvature = function(e, i) -1,
                   psi = list(omega = 0.1, kappa1 = 0.3, pacf1 = 0.8))

# The level tracks the change in the mean.
round(out$eta[c(1, 15, 30, 32, 45, 60)], 3)
#> [1]  0.500 -0.468 -0.181  1.475  2.268  2.079

# The Jacobian that the filter propagates is exact; a central difference
# of the filter agrees with it to the accuracy of the step.
f <- function(v) sum(term_filter(b, rep(0, 60), dd$y,
                                 function(e, i) dd$y[i] - e,
                                 function(e, i) -1,
                                 list(omega = v[1], kappa1 = v[2],
                                      pacf1 = v[3]))$eta)
v0 <- c(0.1, 0.3, 0.8)
h  <- 1e-5
fd <- vapply(1:3, function(j) {
  e <- numeric(3); e[j] <- h
  (f(v0 + e) - f(v0 - e)) / (2 * h)
}, numeric(1))
max(abs(colSums(out$jacobian) - fd))
#> [1] 6.187901e-08

# A panel: the filter runs on each group from its own starting level.
dd$id <- rep(1:3, each = 20)
term_build(gas(p = 1, q = 1, by = id, time = t), dd)
#> <GasTerm> 'gas': score-driven, p = 1, q = 1; 3 group(s)
#>   parameters: omega, kappa1, pacf1

# A developed parameter expands in place, and brings its sub-term's
# penalty with it.
dd$g <- factor(rep(c("u", "v"), 30))
gb <- term_build(gas(p = 1, q = 1, omega ~ ridge(~ g), time = t), dd)
term_params(gb)
#> [1] "omega.(Intercept)" "omega.ridge.gu"    "omega.ridge.gv"   
#> [4] "kappa1"            "pacf1"            
vapply(term_penalties(gb), function(e) e$name, character(1))
#> [1] "omega::ridge(~g)"


# The term is fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  # The series is drawn from the recursion that the term implements,
  # with a gaussian mean.
  set.seed(2)
  om <- 0.05; al <- 0.3; be <- 0.9
  lev <- om / (1 - be)
  yy <- numeric(600)
  for (i in seq_len(600)) {
    yy[i] <- lev + rnorm(1)
    lev <- om + al * (yy[i] - lev) + be * lev
  }
  fd <- data.frame(t = seq_len(600), y = yy)
  ft <- statmodels7::statmod(y ~ 0 + gas(p = 1, q = 1, time = t),
                             distributions7::gaussian1_distrib(), fd)
  # The truth is omega 0.05, kappa 0.3, xi 0.9. xi is reported as the
  # autoregressive coefficient, whereas the estimated coordinate is a
  # partial autocorrelation.
  print(round(coef(ft)$mu, 3))

  # With an intercept the level is held at zero and the intercept
  # estimates omega / (1 - xi) instead: the same model, reparametrized.
  fi <- statmodels7::statmod(y ~ 1 + gas(p = 1, q = 1, time = t),
                             distributions7::gaussian1_distrib(), fd)
  print(c(logLik(ft), logLik(fi)))
  cf <- coef(ft)$mu
  c(intercept  = coef(fi)$mu[["(Intercept)"]],
    from_omega = cf[["gas.omega"]] / (1 - cf[["gas.xi1"]]))
}
#>  gas.omega gas.kappa1    gas.xi1 
#>      0.071      0.329      0.915 
#> [1] -867.599 -867.599
#>  intercept from_omega 
#>  0.8376033  0.8376033 
```
