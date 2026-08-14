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
  distrib = NULL,
  correlated = TRUE,
  label = "random",
  hyper = NULL,
  ...
)
```

## Arguments

- formula:

  A bar formula, `~ lhs | g`, with `g` evaluating to the grouping
  variable in the data.

- distrib:

  The distribution of the effects: `NULL` (the default Gaussian), a
  distributions7 object, or a list of them with one per within-group
  column.

- correlated:

  Logical; whether the default Gaussian lets the within-group effects
  correlate. It is an error together with `distrib`, which says the same
  thing and more.

- label:

  A single non-empty string prefixed to the coefficient names.

- hyper:

  The hyperparameters of the effects' distribution to HOLD, as a named
  vector or list; those not named are estimated. The names are the
  distribution's own parameters, with the within-group column appended
  where there is one copy per column. A name the penalty does not carry
  is reported when the term is built, which is where the penalty first
  exists.

- ...:

  Unused; a named argument here is reported rather than ignored.

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

Two things are said and nothing else is: the formula, and the
distribution of the effects. Which chart the hyperparameters ride, what
they are called, how many there are and where the log-density has a kink
are properties of that distribution, read off it rather than restated
here.

## The distribution of the effects

`distrib` is `NULL`, a distributions7 object, or a list of them with one
per within-group column.

A MULTIVARIATE distribution of the within-group dimension lets the
effects of one group depend on each other, its matrix parameter carrying
the dependence: `mvgaussian_distrib(2, omega = ar1(2))` is a prior whose
precision is autoregressive, `mvstudent_t_distrib(2)` a heavy-tailed
one. Correlation is available exactly for the families that carry a
matrix parameter – a location block as long as the dimension, and a
covariance, precision or scale matrix – which is a property the term
reads rather than a list of admitted names, so a family added later is
covered.

A UNIVARIATE distribution makes the effects independent, the penalty
being the product of the densities. With more than one within-group
column it is a TEMPLATE: one copy per column, each with its own
hyperparameters, since an intercept and a slope are quantities of
different units and a shared scale would price them against each other.
A list of distributions gives one per column explicitly, when the
columns want different priors.

The default is Gaussian: `gaussian1_distrib` at one column, so the
hyperparameter IS the standard deviation of the effects; the
multivariate Gaussian on an unstructured covariance when there are
several and `correlated = TRUE`; the template of the first when
`correlated = FALSE`, one standard deviation per column.

Whatever it is, the distribution is CENTERED, its location parameters
held with
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.html).
A free mean in the effects is confounded with the intercept of the
equation the term sits in, so it is rejected rather than fitted along a
flat direction. The value it is held at is usually zero and is not
policed: it is identified whatever it is, and where the prior is a
transformation of another family the parameter is the mean on the
ORIGINAL scale –
`fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)` is a
log-gamma prior whose own mean is \\\psi(a) - \log a\\, within
\\\sigma^2/2\\ of zero and exactly zero in the limit.

A distribution used as a penalty gives joint-mode (penalized likelihood)
estimation of the effects. With Gaussian effects and a Gaussian response
the Laplace approximation behind a marginal criterion is exact, so the
variance component it returns is the marginal estimate; with any other
prior it is an approximation, and the marginal likelihood is not
computed here.

## The hyperparameters

They are the distribution's own free parameters, and every one of them
is estimated unless it is held. There are two ways to hold one, and they
differ in what is reported rather than in the fit. Holding it inside the
distribution, `fixed(pseudohuber_distrib(), mu = 0, nu = 2)`, removes
it: it becomes a constant of the prior and appears nowhere among the
model's hyperparameters. Naming it in `hyper` keeps it, reported as held
at the value given, which is what a penalized term's own hyperparameter
argument does.

A smooth prior's hyperparameters are estimated by a marginal criterion.
A prior whose log-density has a kink – a Laplace, an elastic net – has
none a marginal criterion can reach, and its hyperparameter is chosen by
a path on a prediction criterion instead.

Every estimated hyperparameter is reported with a standard error and an
interval, shape parameters included. Where one is absent the cause is
the POINT and not a missing derivative: a run that ended where the
criterion has no maximum – a shape escaping towards a limit is the
common case – leaves a curvature of the wrong sign, and no interval
follows from it.

Which PARAMETRIZATION of a family is used matters here in a way it does
not elsewhere. The centred skew normal
([`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.html))
carries the skewness itself, and its map to the direct parametrization
is not twice differentiable at zero skewness: the first derivatives have
a finite limit there and the second ones grow like \\\gamma_1^{-2/3}\\.
A marginal criterion reads the second, and the symmetric bounds put the
starting value at exactly that point, so the direct parametrization
([`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.html))
is the one to use as a prior; its derivatives at \\\alpha = 0\\ are
ordinary numbers.

How well a shape parameter is estimated depends on how many groups there
are, since it is read off that many latent values, and the prior shrinks
them. Measured on effects drawn from a Student t with four degrees of
freedom, twelve observations per group, the prior being a Student t with
\\\nu\\ free: \\\hat\nu\\ is 17.1 at 20 groups, 3.95 at 100 and 4.06 at
500. At 100 the profile has an interior maximum, the criterion falling
from -1618.7 at \\\nu = 4\\ to -1619.3 either side and -1622.0 in the
Gaussian limit. So a shape is worth estimating from a hundred groups or
so and worth holding below that, where it escapes towards the Gaussian
limit and only the scale is really being fitted. A pseudo-Huber's
\\\nu\\ is the weaker case: it is where the loss stops being quadratic
rather than a tail index, and at 40 groups it escapes.

Prediction maps new data onto the levels seen at build time; a level the
term has not seen is rejected.

A random effect is not standardized, and there is no `standardize`
argument to ask for it with; passing one is an error. Its columns are
grouping indicators rather than measured covariates, and its
hyperparameter is a variance component with a meaning of its own.
Dividing each coefficient by the spread of its indicator would weight
the effects by the sizes of the groups, which changes the model rather
than the scale its hyperparameter is read on.

## The block and its penalty

With \\m\\ levels and a within-group design \\Z_i\\ of \\d\\ columns,
the block is the interaction of that design with the group indicators,
ordered group by group,

\$\$Z = \operatorname{diag}(Z_1, \dots, Z_m), \qquad b = (b_1', \dots,
b_m')',\$\$

so the \\d\\ coefficients of one group occupy adjacent positions and the
penalty reads them one block at a time,

\$\$\rho(b; \theta) = -\sum\_{i=1}^m \log f(b_i; \theta),\$\$

for the effects' density \\f\\. Under the default Gaussian that is

\$\$\rho(b; \theta) = \tfrac{1}{2}\sum_i b_i'\Sigma(\theta)^{-1}b_i +
\tfrac{m}{2}\log\lvert \Sigma(\theta)\rvert +
\tfrac{md}{2}\log(2\pi),\$\$

and minimizing the penalized least squares in \\(\beta, b)\\ is the
mixed-model equation at the variance ratio \\\Sigma\\ encodes, so the
minimizer is the best linear unbiased predictor. At \\d = 1\\ and
\\\Sigma = \sigma_b^2\\ it is the ridge, up to the constant that makes
\\\sigma_b\\ estimable.

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

# one variance component, reported as a standard deviation
term_penalty(term_build(random(~ 1 | g), dd))@params
#> [1] "sigma"

# intercepts and slopes, correlated: the covariance of the effects
built <- term_build(random(~ x | g), dd)
term_coef_names(built)
#> [1] "random.a.(Intercept)" "random.a.x"           "random.b.(Intercept)"
#> [4] "random.b.x"           "random.c.(Intercept)" "random.c.x"          
term_penalty(built)@params
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  

# independent, one standard deviation per column
vapply(term_penalties(term_build(random(~ x | g, correlated = FALSE), dd)),
       function(e) e$name, "")
#> [1] "(Intercept)" "x"          

# a heavy-tailed prior, held at four degrees of freedom
t4 <- distributions7::fixed(distributions7::student_t1_distrib(),
                            mu = 0, nu = 4)
term_penalty(term_build(random(~ 1 | g, distrib = t4), dd))@params
#> [1] "sigma"
```
