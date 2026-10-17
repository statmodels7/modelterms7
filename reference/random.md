# Grouped Random-Effect Term

Random intercepts and slopes for a grouping factor: `random(~ 1 | g)`
builds one coefficient per level of `g`, and `random(~ x | g)` one
intercept and one slope per level. The distribution of the effects is
attached as the penalty on those coefficients, which is how a random
effect enters a penalized likelihood.

## Usage

``` r
random(
  formula,
  distrib = NULL,
  correlated = TRUE,
  label = "random",
  hyper = NULL,
  id = NULL,
  ...
)
```

## Arguments

- formula:

  A bar formula, `~ lhs | g`, where `g` evaluates to the grouping
  variable in the data, or `~ lhs | tag | g` with a covariance label in
  the middle. `g` must not be continuous: a double carrying a value that
  is not a whole number is rejected, since its levels would be its
  formatted values.

- distrib:

  The distribution of the effects: `NULL` (the default Gaussian), a
  distributions7 object, or a list of them with one per within-group
  column.

- correlated:

  Logical; whether the default Gaussian lets the within-group effects
  correlate. Giving it together with `distrib` is an error, since
  `distrib` determines the same thing.

- label:

  A single non-empty string prefixed to the coefficient names.

- hyper:

  The hyperparameters of the effects' distribution to hold, as a named
  vector or list; those not named are estimated. For a distribution
  whose log-density has a kink, several values for one name are the grid
  that the path visits, and for any other distribution they are
  rejected. The names are the distribution's own parameters, and the
  within-group column is appended where there is one copy per column. A
  name that the penalty does not carry is rejected when the term is
  built, which is where the penalty first exists.

- id:

  Labels under which the effects' hyperparameters are shared with the
  hyperparameters of other terms carrying the same labels, so that two
  random effects may be given one variance component. The names are
  `hyper`'s, checked at the same point and for the same reason; a single
  unnamed string is accepted where the effects' distribution carries one
  hyperparameter. `NULL`, the default, shares nothing. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md).

- ...:

  Not used. A named argument is rejected with an error naming it, and
  the removed arguments `precision` and `kinks` are rejected with a
  message giving the replacement. An unnamed argument beyond those of
  the signature is rejected as well.

## Value

An object of class
[`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md)
(a specification; see
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The left side of the bar is an ordinary one-sided formula for the
within-group design, with the usual intercept convention: `~ x | g`
carries an intercept and a slope per group and `~ 0 + x | g` carries
only the slope. The block interacts that design with the group
indicators, ordered group by group, so the coefficients of one group are
adjacent.

The formula and the distribution of the effects determine the term. The
parametrization in which the hyperparameters are estimated, their names,
their number and the points where the log-density has a kink are all
properties of that distribution, read from it at build time.

## The covariance label, and the second bar

A formula may carry two bars, `~ 1 + x | u | g`, following brms. The
last position is the grouping variable and the middle one is a label:
terms carrying the same label and the same grouping describe effects
that are correlated with each other, so an intercept in one equation and
a slope in another can share one covariance block. R nests bars to the
left, so the parser itself yields that reading.

The label is a name and not data. It is read as a symbol and never
evaluated, and a column of the data that carries the same name is not
read, so `~ 1 | u | g` builds the same block as `~ 1 | g`, whatever `u`
is.

A labelled term declares no penalty of its own. The covariance block
spans columns to which the term has no access, possibly in another
equation, so the prior over it belongs to the class and is built by the
fitting layer.
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
returns an empty list and
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
is `NULL`. Consequently `hyper` and `id` are rejected when the term is
built, because the hyperparameters are those of the class, and a
`distrib` that is given on a labelled term is the class's joint prior:
it must be multivariate, and its dimension is checked where the class is
known.

The two arguments have different scopes. `correlated` indicates whether
the columns of this term may depend on each other, and the label
indicates that its columns depend on those of other terms. The
combination of `correlated = FALSE` with a label is contradictory and is
rejected when the term is built.

## The distribution of the effects

`distrib` is `NULL`, a distributions7 object, or a list of them with one
per within-group column.

A multivariate distribution of the within-group dimension lets the
effects of one group depend on each other, and its matrix parameter
carries the dependence: `mvgaussian2_distrib(2, ar1(2))` is a prior
whose precision is autoregressive, `mvstudent_t1_distrib(2)` a
heavy-tailed one. Correlated effects require a family that carries a
matrix parameter, that is, a location block as long as the dimension
together with a covariance, precision or scale matrix, and that defines
the mixed response-parameter derivative block read by a marginal
criterion. The term reads these properties from the family, so a family
added later is covered without an edit here.

A univariate distribution makes the effects independent, and the penalty
is the product of the densities. With more than one within-group column
it is a template: one copy per column, each with its own
hyperparameters, since an intercept and a slope are quantities of
different units and a shared scale would price them against each other.
A list of distributions gives one per column explicitly, when the
columns need different priors.

The default is Gaussian: `gaussian1_distrib` at one column, so the
hyperparameter is the standard deviation of the effects; the
multivariate Gaussian on an unstructured covariance when there are
several columns and `correlated = TRUE`; and the template of the first
when `correlated = FALSE`, with one standard deviation per column. The
default has its location parameters held at zero.

A distribution supplied in `distrib` must have its location parameters
held with
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html).
A free location in the effects is confounded with the intercept of the
equation in which the term sits, so a free location is rejected at build
time with a message naming it. The value at which it is held is usually
zero and is not checked, since the model is identified whatever the
value is. Where the prior is a transformation of another family, the
parameter is the mean on the original scale, so
`fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)` is a
log-gamma prior with shape \\a = 1/\sigma^2\\, where \\\sigma^2\\ is the
parameter `sigma2` of the gamma. The mean of that prior, \\\psi(a) -
\log a\\, lies between \\-\sigma^2\\ and \\-\sigma^2/2\\ and tends to
zero as \\\sigma^2\\ tends to zero.

A distribution used as a penalty gives joint-mode (penalized likelihood)
estimation of the effects. With Gaussian effects and a Gaussian response
the Laplace approximation behind a marginal criterion is exact, so the
variance component that the criterion returns is the marginal estimate.
With any other prior the Laplace approximation is not exact, and
modelterms7 does not compute the marginal likelihood.

## The hyperparameters

They are the distribution's own free parameters, and every one of them
is estimated unless it is held. There are two ways to hold one, and the
fit is the same either way, but the reported output differs. Holding it
inside the distribution, `fixed(pseudohuber_distrib(), mu = 0, nu = 2)`,
removes it: it becomes a constant of the prior and is not among the
model's hyperparameters. Naming it in `hyper` keeps it, reported as held
at the value given, as a penalized term's own hyperparameter argument
does.

The hyperparameters of a prior with a smooth log-density are estimated
by the outer criterion. A prior whose log-density has a kink, a Laplace
or an elastic net, has none that a marginal criterion can reach, and its
hyperparameter is chosen by a path over the sparse criterion of
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html),
[`statmodels7::bic()`](https://statmodels7.github.io/statmodels7/reference/aic.html)
by default.

An estimated hyperparameter, including a shape parameter, is reported
with a standard error and an interval. The standard error is absent
where the criterion at the point at which the run ended has a curvature
of the wrong sign, since such a point does not yield an interval. A
shape parameter that moves toward a limit is the common case.

The parametrization of a family that is used as a prior matters. The
centered skew normal
([`distributions7::skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.html))
carries the skewness itself, and its map to the direct parametrization
is not twice differentiable at zero skewness. The derivatives of the
log-density with respect to the skewness have a finite first-order limit
there, and the second-order ones grow like \\\gamma_1^{-2/3}\\. A
marginal criterion reads the second-order ones, and the symmetric bounds
place the starting value at that point, so the direct parametrization
([`distributions7::skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.html))
is the one to use as a prior; its derivatives at \\\alpha = 0\\ are
finite.

How well a shape parameter is estimated depends on the number of groups,
since the parameter is read from that many latent values, which the
prior shrinks. With few groups the estimate of a shape parameter can
move to the Gaussian limit, and then only the scale is estimated. With
more groups the estimate stays finite, but the criterion is flat over
the small values, so the estimate locates the shape only to its order of
magnitude. A shape parameter can be held with
[`distributions7::fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.html)
when the number of groups is small, and an estimate from a moderate
number of groups describes the weight of the tail without fixing its
value. The \\\nu\\ of a pseudo-Huber prior, the point at which the loss
stops being quadratic, is subject to the same limitation.

Prediction maps new data onto the levels seen at build time. A level
that was not present at build time, or a missing grouping value, is
rejected by
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
unless `unseen = "zero"` is given, in which case those rows have no
effect of their own.

A random effect is not standardized, and `random()` has no `standardize`
argument; passing one is an error. Its columns are grouping indicators
and its hyperparameter is a variance component with a meaning of its
own. Dividing each coefficient by the spread of its indicator would
weight the effects by the sizes of the groups, which changes the model
itself.

## The block and its penalty

With \\m\\ levels and a within-group design \\Z_i\\ of \\d\\ columns,
the block is the interaction of that design with the group indicators,
ordered group by group. When the rows of the data are ordered by group,
the block and the coefficient vector are

\$\$Z = \operatorname{diag}(Z_1, \dots, Z_m), \qquad b = (b_1', \dots,
b_m')',\$\$

and otherwise the rows of this matrix are permuted to follow the data.
The \\d\\ coefficients of one group occupy adjacent positions and the
penalty reads them one block at a time,

\$\$\rho(b; \theta) = -\sum\_{i=1}^m \log f(b_i; \theta),\$\$

for the effects' density \\f\\. Under the default Gaussian that is

\$\$\rho(b; \theta) = \tfrac{1}{2}\sum_i b_i'\Sigma(\theta)^{-1}b_i +
\tfrac{m}{2}\log\lvert \Sigma(\theta)\rvert +
\tfrac{md}{2}\log(2\pi),\$\$

and minimizing the penalized least squares in \\(\beta, b)\\ is the
mixed-model equation at the variance ratio that \\\Sigma\\ encodes, so
the minimizer is the best linear unbiased predictor. At \\d = 1\\ and
\\\Sigma = \sigma_b^2\\ it is the ridge, up to the constant that makes
\\\sigma_b\\ estimable.

## Smooth random effects

Where the left side of the bar is an
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) or a
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)
call, `random(~ s(x) | g)`, each level of `g` has a random curve of its
own (a random surface for
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)),
and the term is a
[`RandomSmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomSmoothTerm.md).
When the smooth's columns are split into the null space of its penalty
\\N(x)\\ and the penalized coordinates \\Z(x)\\, group \\i\\ contributes

\$\$g_i(x) = a\_{0i} + N(x)a_i + Z(x)\delta_i, \qquad (a\_{0i}, a_i)
\sim \mathrm{N}(0, \Sigma), \qquad \delta_i \sim \mathrm{N}(0, \sigma^2
I),\$\$

independently over the groups. The null space is whatever the smoother's
operator leaves unpenalized, less the constant: a straight line for a
second-derivative penalty, a sine and a cosine for the harmonic
operator, nothing for a cyclic smoother. The argument `null_space` of
the smoother determines what the null space contributes: under `"keep"`
its columns carry \\\Sigma\\ with the random level, under `"shrink"`
they join the penalized coordinates, and under `"drop"` they are absent.
`~ 0 + s(x) | g` removes the random level.

\\\Sigma\\ is the prior that `distrib` and `correlated` describe for the
within-group columns of an ordinary random effect, here over the random
level and the null space. With an unstructured \\\Sigma\\ the model does
not depend on the origin and the scale of `x`, nor on the phase of a
periodic one. The penalized coordinates are whitened, so that their
prior is independent with one standard deviation, reported as `sigma` of
the entry `smooth`. A smoother with several smoothing parameters (an
adaptive smoother, an anisotropic tensor product) keeps them: the
precision is \\\sum_k \lambda_k S_k\\ in every group, and the summary
reports \\\lambda_k^{-1/2}\\. All the groups share these
hyperparameters, so the curves have one roughness. The population curve
is a separate
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) term
in the same formula, as in `y ~ s(x) + random(~ s(x) | g)`.

The smooth is written without `by`, `sparse`, `hyper` or `id`. The
hyperparameters belong to the random effect, and `hyper` and `id` of
this function name them (`sigma.smooth`, `lambda1.smooth`, ...). A
covariance label between two bars is rejected.

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

# The term has one variance component, which is reported as a standard
# deviation.
term_penalty(term_build(random(~ 1 | g), dd))@params
#> [1] "sigma"

# Intercepts and slopes are correlated, and the covariance of the effects
# is the hyperparameter.
built <- term_build(random(~ x | g), dd)
term_coef_names(built)
#> [1] "random.a.(Intercept)" "random.a.x"           "random.b.(Intercept)"
#> [4] "random.b.x"           "random.c.(Intercept)" "random.c.x"          
term_penalty(built)@params
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  

# The effects are independent, and each column has its own standard
# deviation.
vapply(term_penalties(term_build(random(~ x | g, correlated = FALSE), dd)),
       function(e) e$name, "")
#> [1] "(Intercept)" "x"          

# A heavy-tailed prior is held at four degrees of freedom.
t4 <- distributions7::fixed(distributions7::student_t1_distrib(),
                            mu = 0, nu = 4)
term_penalty(term_build(random(~ 1 | g, distrib = t4), dd))@params
#> [1] "sigma"

# A covariance label is read from the middle of two bars and is reported by
# term_tag(). It is a name: the block is the one bar formula's.
term_tag(random(~ 1 | u | g))
#> [1] "u"
identical(term_coef_names(term_build(random(~ 1 | u | g), dd)),
          term_coef_names(term_build(random(~ 1 | g), dd)))
#> [1] TRUE

# A continuous grouping variable is rejected rather than turned into one
# level per row.
try(term_build(random(~ 1 | x), dd))
#> Error : the grouping variable 'x' is continuous, and a grouping variable's
#>   values are labels rather than measurements: its levels would be the
#>   formatted numbers, one group per distinct value. Write factor(x) if
#>   they really are labels.

# A smooth random effect has one random curve per level of g.
sd <- data.frame(x = runif(60), g = factor(rep(c("a", "b", "c"), 20)))
rs <- term_build(random(~ s(x, smoother = basis7::bspline_smooth(k = 5)) |
                          g), sd)
term_group(rs)$names
#> [1] "(Intercept)" "lin"         "z1"          "z2"          "z3"         


# The model is fitted to data simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(6)
  fd <- data.frame(gr = factor(rep(1:20, each = 15)), x = rnorm(300))
  bb <- rnorm(20, sd = 0.8)
  fd$y <- 1 + 0.5 * fd$x + bb[fd$gr] + rnorm(300, sd = 0.5)
  cf <- coef(statmodels7::statmod(y ~ x + random(~1 | gr),
                                  distributions7::gaussian1_distrib(), fd))$mu
  # The true slope is 0.5, and the shrunken effects track the ones drawn.
  round(c(slope = cf[["x"]],
          cor = cor(cf[grep("^random", names(cf))], bb)), 3)
}
#> slope   cor 
#> 0.472 0.989 
```
