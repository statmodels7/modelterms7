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
gas(
  p = 1,
  q = 1,
  by = NULL,
  time = NULL,
  deviations = FALSE,
  penalty = c("none", "lasso", "ridge"),
  label = "gas"
)
```

## Arguments

- p:

  The number of score lags. Defaults to 1.

- q:

  The number of autoregressive lags. Defaults to 1.

- by:

  An optional grouping variable, evaluated in the data; each group is
  filtered independently, from its own starting level.

- time:

  An optional ordering variable, evaluated in the data.

- deviations:

  Whether each group carries a deviation from the population parameters:
  `FALSE` (default), `TRUE` for every parameter, or a character vector
  naming the parameters that carry one. Requires `by`.

- penalty:

  One of `"none"` (default), `"lasso"` or `"ridge"`, applied to the
  deviations. Requires them.

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
\dots, a_p\\, and the persistence. The persistence is carried by
**partial autocorrelations** rather than by the coefficients \\b_j\\:
the stationary region of an autoregression is not a box, so no
collection of scalar links covers it, while the partial autocorrelations
each range over \\(-1, 1)\\ independently and the Levinson-Durbin
recursion carries them onto the coefficients bijectively. At \\q = 1\\
the two coincide. The coordinate is named for the chart it lives on,
`pacf1` and so on, following the convention of parameters7.

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

### A population value and a deviation per group

By default every group of a panel is filtered with the same parameters.
`deviations` gives each group its own, written as a population value and
a departure from it, \$\$\psi\_{j,i} = g_j^{-1}\\\left(g_j(\psi_j) +
\delta\_{j,i}\right),\$\$ the deviation acting on the unconstrained
scale of the chart the parameter lives on, so that a persistence stays
inside \\(-1, 1)\\ whatever the deviation is. The deviations are
parameters of the term, named `omega.dev.1` and so on after the
parameter and the level, and they carry the identity link, being
unconstrained already.

They are parameters and not a penalty on the per-group values through a
difference matrix, which is what the same model looks like written the
other way. The difference decides what can be fitted: a penalty over a
general map is the generalized-lasso problem, whose proximal operator
does not split by coordinate, whereas a deviation named as a coordinate
is reached by a soft threshold and by a coordinate descent unchanged.
`penalty` shrinks them towards zero, which is towards a panel that is
homogeneous in that parameter, and `"lasso"` sets the deviations of the
groups that do not need one exactly to it.

The penalty is also what identifies them. A parameter and its \\m\\
deviations are \\m+1\\ numbers describing \\m\\ group values, so adding
a constant to \\g_j(\psi_j)\\ and subtracting it from every
\\\delta\_{j,i}\\ leaves the filter and its likelihood exactly
unchanged: without a penalty on the deviations the likelihood is flat
along one direction per parameter carrying them. This is the
parametrization of a random effect, and it is identified in the same way
– there by a variance component, here by the penalty, which selects the
deviations of smallest size among those that describe the same panel.
`penalty = "none"` is therefore for reading a filter at given parameters
rather than for fitting one.

The parameters a specification reports are the population ones alone:
how many groups there are is a property of the data, so the deviations
appear once the term is built.

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
#> [1] "omega" "a1"    "pacf1" "pacf2"
```
