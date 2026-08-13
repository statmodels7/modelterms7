# Score-Driven Dynamics

A generalized autoregressive score component (creal2013, harvey2013): a
level \\f_t\\ added to the linear predictor, driven by the score of the
observation density at the previous times, \$\$f_t = \omega +
\sum\_{i=1}^{p} a_i s\_{t-i} + \sum\_{j=1}^{q} b_j f\_{t-j},\$\$ with
\\s_t = \partial \ell_t / \partial \eta_t\\ the derivative of the
log-likelihood contribution with respect to the predictor it is
evaluated at.

## Usage

``` r
gas(p = 1, q = 1, ..., by = NULL, time = NULL, links = NULL, label = "gas")
```

## Arguments

- p:

  The number of score lags. Defaults to 1.

- q:

  The number of autoregressive lags. Defaults to 1.

- ...:

  Two-sided formulas whose left side names a parameter, one per
  parameter to be developed with covariates, e.g. `alpha1 ~ s(x)`; see
  the section above.

- by:

  An optional grouping variable, evaluated in the data; each group is
  filtered independently, from its own starting level. A FORMULA here is
  the shorthand giving the same subformula to every parameter, and then
  no grouping is implied.

- time:

  An optional ordering variable, evaluated in the data.

- links:

  An optional named list of linkfunctions7 links over the parameters of
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
  overriding the defaults described above. A deviation cannot be named:
  it is unconstrained by construction, acting on the scale its
  parameter's own link defines.

- label:

  A single non-empty string naming the term.

## Value

An object of class
[`GasTerm`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The term adds no columns. The predictor at one time depends on the data
at the previous ones, so the contribution cannot be written as a block,
and
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
runs the recursion instead. That is what makes it a
[`structural_term`](https://statmodels7.github.io/modelterms7/reference/structural_term.md).

What drives the recursion is the score of whatever distribution the
model carries, so the same term is a GARCH-like volatility model when it
enters the scale of a Gaussian, a dynamic count model when it enters the
mean of a Poisson, and a robust location filter when it enters a Student
t: a heavy-tailed score is bounded in the observation, so an outlier
moves the level by a bounded amount rather than in proportion to its
size.

### The parameters and their chart

The parameters are the level \\\omega\\, the score loadings \\a_1,
\dots, a_p\\, and the persistence. Each is estimated on the
unconstrained scale of a link, and `links` overrides any of them; the
defaults are the following.

The level carries the identity, being unconstrained. The loadings carry
the **log** link: a positive loading responds in the direction of the
score, which is the case the score-driven literature writes, and
positivity is then structural – a deviation or a submodel moves the
loading on the log scale and no group or observation can take a negative
one. A loading that must be free in sign is asked for with
`links = list(alpha1 = linkfunctions7::identity_link())`.

The persistence is carried by **partial autocorrelations** rather than
by the coefficients \\b_j\\: the stationary region of an autoregression
is not a box, so no collection of scalar links covers it, while the
partial autocorrelations each range over \\(-1, 1)\\ independently and
the Levinson-Durbin recursion carries them onto the coefficients
bijectively. At \\q = 1\\ the two coincide. The coordinate is named for
the chart it lives on, `pacf1` and so on, following the convention of
parameters7.

Whatever the links, a parameter modeled per group or per observation
stays in its own set: a departure acts on the unconstrained scale, so a
loading on the log link is positive and a persistence on the rhobit link
is stationary at every observation, whatever the departure is.

### One parameter at a time

The level this term drives is a scalar, so it enters the predictor of
one distribution parameter. In the general formulation the level is a
vector with one entry per modeled parameter, \\\omega\\ a vector and the
loadings \\A_i\\ and \\B_j\\ matrices, which lets the scale respond to
the score of the location and the other way round. The recursion
generalizes mechanically, and so does the derivative propagated with it;
what is missing is a way to say that one filter spans several
parameters, since a term written inside the formula of one parameter has
no place to declare it. That belongs to the model layer. The persistence
would also need a different chart: the partial-autocorrelation
construction below is a scalar one, and the stationary region of a
matrix autoregression is a bound on the spectral radius of its companion
matrix rather than a box.

The score driving the recursion is used unscaled. The general
formulation carries a scaling matrix, usually an inverse information,
which the curvature this term already receives would supply.

### Groups and time

`by` filters each group independently, which is what a panel of short
series needs, and `time` gives the order within a group. Without `time`
the rows are taken in the order they appear.

### A parameter developed with covariates

A two-sided formula in `...` whose left side names a parameter develops
it as \\\psi\_{j,t} = g_j^{-1}(z_t^\top\gamma_j)\\, the design \\Z\\
built from the right-hand side through
[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
so it takes any additive term of the package:

    gas(p = 1, q = 1, omega ~ ridge(~g), alpha1 ~ s(x),
        pacf1 ~ random(~1 | id), by = id)

The development acts on the unconstrained scale of the parameter's own
link, which is what keeps every per-observation value in the parameter's
own set whatever the coefficients are: a loading on the log link is
positive at every observation, a persistence on the rhobit chart is
inside \\(-1, 1)\\ at every observation, and at \\q = 1\\ that bounds
the recursion's growth step by step. The coefficients \\\gamma_j\\ are
the term's parameters, unconstrained and on the identity link; the
penalties the sub-terms carry are reported through
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
under the key `parameter::subterm`. A parameter that varies by
observation changes the recursion itself, \$\$f_t = \omega_t + \sum_i
a\_{i,t}\\ s\_{t-i} + \sum_j b\_{j,t}\\ f\_{t-j},\$\$ with \\b_t\\ from
the Levinson-Durbin map of that observation's partial autocorrelations,
and the filter, its derivative, the reverse pass and the curvature all
run the general recursion.

`by = ~f` (a formula, where a grouping variable is a bare symbol) is the
shorthand giving the same subformula to every parameter; mixing it with
per-parameter formulas is an error. A structural term, and a term whose
block moves with its own coefficients, are rejected inside a subformula.

### A population value and a departure per group

The panel case is one subformula:
`gas(omega ~ random(~1 | id), by = id)` is a population value (the
intercept of the development) plus one unconstrained departure per
group, shrunk by the random intercept's own ridge, whose hyperparameter
a fitting layer estimates. `lasso(...)` in the subformula sets the
departures of the groups that do not need one exactly to zero.

The departures act on the unconstrained scale of the parameter's chart,
so a persistence stays inside \\(-1, 1)\\ whatever the departure is.
Their penalty is also what identifies them: a population value and \\m\\
departures are \\m+1\\ numbers describing \\m\\ group values, so without
the penalty the likelihood is flat along one direction per developed
parameter. This is the parametrization of a random effect and it is
identified the same way, there by a variance component, here by the
penalty. An unpenalized development over group indicators is therefore
for reading a filter at given parameters rather than for fitting one.

Earlier releases spelled this case as `deviations =` with a `penalty =`;
both arguments are gone, the subformula reproducing them exactly (the
same fit to the printed digit, hyperparameter included) while covering
what they could not.

## References

Creal, D., Koopman, S. J. and Lucas, A. (2013). Generalized
autoregressive score models with applications. *Journal of Applied
Econometrics*, 28(5), 777–795.

Harvey, A. C. (2013). *Dynamic Models for Volatility and Heavy Tails*.
Cambridge University Press.

## See also

[`regime`](https://statmodels7.github.io/modelterms7/reference/regime.md)

## Examples

``` r
term_params(gas(p = 1, q = 2))
#> [1] "omega"  "alpha1" "pacf1"  "pacf2" 
```
