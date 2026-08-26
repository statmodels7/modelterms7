# Grouped Random-Effect Term

Random intercepts and slopes for a grouping factor: `random(~ 1 | g)`
builds one coefficient per level of `g`, and `random(~ x | g)` one
intercept and one slope per level, with the distribution of the effects
attached as the penalty on those coefficients. That is what a random
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

  The hyperparameters of the effects' distribution to hold, as a named
  vector or list; those not named are estimated. The names are the
  distribution's own parameters, with the within-group column appended
  where there is one copy per column. A name the penalty does not carry
  is reported when the term is built, which is where the penalty first
  exists.

- ...:

  Unused. A named argument here is reported by name, so a removed one
  such as `precision` or `kinks` gets a message saying what replaced it.

## Value

An object of class
[`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md)
(a specification; see
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The left side of the bar is an ordinary one-sided formula for the
within-group design, with the usual intercept convention: `~ x | g`
carries an intercept and a slope per group and `~ 0 + x | g` the slope
alone. The block interacts that design with the group indicators,
ordered group by group, so the coefficients of one group are adjacent.

The constructor asks for two things: the formula and the distribution of
the effects. Which chart the hyperparameters ride, what they are called,
how many there are and where the log-density has a kink are all
properties of that distribution, read off it at build time.

## The distribution of the effects

`distrib` is `NULL`, a distributions7 object, or a list of them with one
per within-group column.

A multivariate distribution of the within-group dimension lets the
effects of one group depend on each other, its matrix parameter carrying
the dependence: `mvgaussian_distrib(2, omega = ar1(2))` is a prior whose
precision is autoregressive, `mvstudent_t_distrib(2)` a heavy-tailed
one. Correlation is available exactly for the families that carry a
matrix parameter: a location block as long as the dimension, together
with a covariance, precision or scale matrix. The term reads that
property off the family, so a family added later is covered without an
edit here.

A univariate distribution makes the effects independent, the penalty
being the product of the densities. With more than one within-group
column it is a template: one copy per column, each with its own
hyperparameters, since an intercept and a slope are quantities of
different units and a shared scale would price them against each other.
A list of distributions gives one per column explicitly, when the
columns want different priors.

The default is Gaussian: `gaussian1_distrib` at one column, so the
hyperparameter IS the standard deviation of the effects; the
multivariate Gaussian on an unstructured covariance when there are
several and `correlated = TRUE`; the template of the first when
`correlated = FALSE`, one standard deviation per column.

Whatever it is, the distribution is centered, its location parameters
held with
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html).
A free mean in the effects is confounded with the intercept of the
equation the term sits in, so a free location is rejected at build time
with a message naming it. The value it is held at is usually zero and is
not policed: the model is identified whatever it is. Where the prior is
a transformation of another family the parameter is the mean on the
original scale, so
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
is estimated unless it is held. There are two ways to hold one, and the
fit is the same either way; what differs is what gets reported. Holding
it inside the distribution,
`fixed(pseudohuber_distrib(), mu = 0, nu = 2)`, removes it: it becomes a
constant of the prior and appears nowhere among the model's
hyperparameters. Naming it in `hyper` keeps it, reported as held at the
value given, as a penalized term's own hyperparameter argument does.

A smooth prior's hyperparameters are estimated by a marginal criterion.
A prior whose log-density has a kink, a Laplace or an elastic net, has
none a marginal criterion can reach, and its hyperparameter is chosen by
a path on a prediction criterion instead.

Every estimated hyperparameter is reported with a standard error and an
interval, shape parameters included. Where one is absent the cause is
the point the run ended at: a criterion with no maximum there leaves a
curvature of the wrong sign, and no interval follows from it. A shape
escaping toward a limit is the common case.

Which parametrization of a family is used matters here in a way it does
not elsewhere. The centered skew normal
([`distributions7::skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.html))
carries the skewness itself, and its map to the direct parametrization
is not twice differentiable at zero skewness: the first derivatives have
a finite limit there and the second ones grow like \\\gamma_1^{-2/3}\\.
A marginal criterion reads the second, and the symmetric bounds put the
starting value at exactly that point, so the direct parametrization
([`distributions7::skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.html))
is the one to use as a prior; its derivatives at \\\alpha = 0\\ are
ordinary numbers.

How well a shape parameter is estimated depends on how many groups there
are, since it is read off that many latent values and the prior shrinks
them. Measured on effects drawn from a standard Student t with four
degrees of freedom, twelve observations per group and unit residual
standard deviation, the prior being a Student t with \\\nu\\ free and
the criterion
[`statmodels7::reml()`](https://statmodels7.github.io/statmodels7/reference/reml.html):

|        |             |                |
|--------|-------------|----------------|
| groups | \\\hat\nu\\ | \\\hat\sigma\\ |
| 20     | 5.97e+04    | 0.769          |
| 100    | 1.98        | 0.817          |
| 500    | 2.65        | 0.923          |

At twenty groups the shape escapes to the Gaussian limit and only the
scale is really being fitted. From a hundred it stays finite, and the
profile is decisive about that much: with \\\nu\\ held, the criterion is
-1922.4 at 3, -1923.7 at 4, -1924.9 at 5 and -1936.1 in the Gaussian
limit. What it is not decisive about is the value, the profile being
flat enough over the small integers that a single sample locates \\\nu\\
to little better than its order of magnitude. Estimate a shape from a
hundred groups or so, hold it below that, and read the estimate as a
statement about the tail rather than a measurement of it.

A pseudo-Huber's \\\nu\\ is the weaker case, being the point at which
the loss stops being quadratic; at 40 groups it escapes.

Prediction maps new data onto the levels seen at build time; a level the
term has not seen is rejected.

A random effect is not standardized, and there is no `standardize`
argument to ask for it with; passing one is an error. Its columns are
grouping indicators and its hyperparameter is a variance component with
a meaning of its own. Dividing each coefficient by the spread of its
indicator would weight the effects by the sizes of the groups, which
changes the model itself.

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

[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md),
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)

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


# Fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(6)
  fd <- data.frame(gr = factor(rep(1:20, each = 15)), x = rnorm(300))
  bb <- rnorm(20, sd = 0.8)
  fd$y <- 1 + 0.5 * fd$x + bb[fd$gr] + rnorm(300, sd = 0.5)
  cf <- coef(statmodels7::statmod(y ~ x + random(~1 | gr),
                                  distributions7::gaussian1_distrib(), fd))$mu
  # truth: a slope of 0.5, and the shrunken effects track the ones drawn
  round(c(slope = cf[["x"]],
          cor = cor(cf[grep("^random", names(cf))], bb)), 3)
}
#> slope   cor 
#> 0.471 0.989 
```
