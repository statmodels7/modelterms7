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
gas(p = 1, q = 1, by = NULL, time = NULL, label = "gas")
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

### Groups and time

`by` filters each group independently, which is what a panel of short
series needs, and `time` gives the order within a group. Without `time`
the rows are taken in the order they appear.

## References

Creal, D., Koopman, S. J. and Lucas, A. (2013). Generalized
autoregressive score models with applications. *Journal of Applied
Econometrics*, 28(5), 777–795.

Harvey, A. C. (2013). *Dynamic Models for Volatility and Heavy Tails*.
Cambridge University Press.

## Examples

``` r
term_params(gas(p = 1, q = 2))
#> [1] "omega" "a1"    "pacf1" "pacf2"
```
