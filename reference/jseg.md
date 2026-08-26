# Segmented-with-Jump Term: Slope and Level Both Changing

A covariate whose slope *and* level both change at \\K\\ break-points
estimated with everything else, the union of
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
at the same points. Written in the equation of any parameter of any
distribution, the term contributes

\$\$\beta x_i + \sum\_{k=1}^{K} \bigl\[\delta_k\\\mathbb{1}(x_i \geq
\psi_k) + \gamma_k\\(x_i - \psi_k)\_{+}\bigr\]\$\$

to that parameter's linear predictor, so that in a model \\g(\theta_i) =
\eta_i\\ with the rest of the equation supplying \\z_i'\alpha\\,

\$\$\eta_i = z_i'\alpha + \beta x_i + \sum\_{k=1}^{K}
\bigl\[\delta_k\\\mathbb{1}(x_i \geq \psi_k) + \gamma_k\\(x_i -
\psi_k)\_{+}\bigr\].\$\$

\\\beta\\ is the slope before the first break-point, \\\gamma_k\\ the
change of slope at \\\psi_k\\ and \\\delta_k\\ the change of level
there.

## Usage

``` r
jseg(
  x,
  ...,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  c0 = 0.05,
  smoothed = NULL,
  marginal = FALSE,
  n_boot = 10,
  label = "jseg"
)
```

## Arguments

- x:

  The covariate, an expression evaluated in the data.

- ...:

  Two-sided formulas developing the term's own coefficients; see the
  section above. Cannot be combined with `by`.

- npsi:

  The number of break-points. Defaults to 1.

- psi:

  Optional starting positions; defaults to evenly spaced quantiles of
  the covariate. Where a break-point carries a development they seed it,
  each starting vector solving \\Wp_k \approx \psi_k^{0}\\.

- by:

  A one-sided formula giving every coefficient of the term the same
  development, e.g. `by = ~0 + g` for an independent set per level of a
  factor. A bare variable is rejected; write the formula.

- linear:

  Whether the block carries the linear effect \\\beta x\\. Defaults to
  `TRUE`.

- c0:

  The starting value of the scaling factor that separates the
  observations from the break-point, as a fraction of the distance to
  the ends of the range. Defaults to `0.05`, the value fasola2018
  recommend; see Details.

- smoothed:

  `NULL` (the default: the construction exactly as documented above) or
  a penalties7
  [`penalties7::abs_smoother()`](https://statmodels7.github.io/penalties7/reference/abs_smoother.html),
  e.g.
  [`penalties7::smooth_probit()`](https://statmodels7.github.io/penalties7/reference/smooth_probit.html).
  The smoother replaces the step and the hinge by their smooth versions,
  \\(1 + s'(u))/2\\ and \\(u + s(u))/2\\, so every break-point becomes
  an ordinary parameter of a \\C^\infty\\ model: there is no working
  parametrization, no auxiliary coefficient and no scaling schedule
  (`c0` is ignored, with a message), the block is the true Jacobian and
  the term is fitted by Gauss-Newton like
  [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md). A
  development of a break-point, `psi ~ random(~1 | id)` and penalized
  ones included, is then legal for every kind, the read-off that
  constrained the discontinuous constructions having gone. The
  smoother's width is resolved at build from the covariate's spacing
  (the median gap between distinct values, within groups where a
  break-point development supplies a partition) unless the object
  carries one, and is reported: it is the width of the transition, the
  bent-cable reading, and the smoothing bias it buys is confined to a
  window of that width (probit, quintic) or decays as \\c/(4\|u\|)\\
  (hyperbolic). The objective is still multimodal in the positions –
  smoothing rounds the local optima and does not remove them, so the
  profile start and the `n_boot` restarts stay necessary; a smoothed fit
  from a bad start has been measured converging to an absurd local
  optimum while reporting success.

- marginal:

  Whether the break-point is a latent variable per group, integrated out
  of the likelihood exactly instead of being an estimated position.
  `FALSE`, the default, is the construction documented above. `TRUE`
  requires the subformula `psi ~ random(~1 | g)` and returns a
  structural term of the likelihood shape
  ([`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md));
  see the section of
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
  whose step model is where the marginal buys the most: for a `seg` term
  the native random-changepoint fit (`psi ~ random(~1 | g)` without
  `marginal`) is measured equivalent and remains the recommended route,
  the marginal being the exact-likelihood alternative. One break-point
  here: the conditional is smooth in the position, so the integral runs
  on a Gauss-Kronrod panel per interval between a group's ordered
  observations, and several latents would need a product quadrature
  whose component count no fitting layer can carry.

- n_boot:

  How many bootstrap restarts the fitting layer runs after the iteration
  first converges (Wood 2001, the device `segmented` runs by default):
  each restart re-estimates on a bootstrap resample from the current
  break-points and then on the data again from where the resample ended,
  keeping the better fit. The objective has local optima in the
  break-points, and with several of them a grid start does not reach the
  right basin from every placement. Defaults to 10, `segmented`'s own
  default; 0 disables. The term itself only declares the value; running
  the restarts is the fitting layer's, as with a penalty's
  hyperparameters.

- label:

  A single non-empty string prefixed to the coefficient names.

## Value

An object of class
[`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
(a specification; see
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)),
or of class
[`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md)
with `marginal = TRUE`.

## Details

Both devices of the two constructions are used at once: the truncated
line is differentiable in the break-point and the step is linearized by
the identity of fasola2018, on the rescaled covariate its scaling
schedule provides (see
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)).
What is particular to the joint term is the read-off, since the
truncated line depends on the break-point as well and reading only the
step's pair would discard what the change of slope says about it.
Linearizing both parts about the previous position, and writing \\h\\
for the increment, \\U\\ for the truncated line and \\I = Z -
\psi^{0}W\\ for the indicator,

\$\$\gamma\\(x-\psi)\_{+} \approx \gamma U - \gamma h I, \qquad
\delta\\\mathbb{1}(x\>\psi) = \delta Z - \delta\psi W,\$\$

so the fitted coefficients of \\Z\\ and \\W\\ are \\a = \delta - \gamma
h\\ and \\b = \gamma h\psi^{0} - \delta\psi^{0} - \delta h\\, and the
increment solves

\$\$\gamma h^{2} + a h + (b + a\psi^{0}) = 0,\$\$

the root of smaller modulus. At \\\gamma = 0\\ the quadratic degenerates
to \\h = -(b + a\psi^{0})/a\\, that is \\\psi = -b/a\\, so the pure step
is a case this contains, never an exception to it.

Everything
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
documents about the scaling schedule, the confinement of a break-point
and the choice of starting positions applies here unchanged.

## The coefficients

`beta` (present when `linear = TRUE`), `gamma1` ... `gammaK`, the
changes of slope, then `delta1` ... `deltaK`, the changes of level, then
`g1` ... `gK`, the auxiliary coefficients from which the break-points
are read. Ask
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for the break-points.

## Developing a coefficient with covariates

`beta` and `gammaK` take a development on any design, both entering the
block linearly. `deltaK` and `psiK` enter the quadratic above, which is
a scalar identity: it splits observation by observation only where the
design has one column per group and each observation belongs to one, so
`jseg(x, by = ~0 + g)` is an independent set of everything per level and
a development of `delta` or `psi` on any other design is rejected. The
componentwise reading that would remain diverges whenever a change of
level passes near zero mid-iteration, which on a joint model it
routinely does, so it is rejected instead of shipped diverging.

## The marginal construction

`jseg(x, psi ~ random(~1 | g), marginal = TRUE)` integrates a latent
break-point per group out of the likelihood exactly, on the
Gauss-Kronrod panels of
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)'s
marginal with the change of level entering each interval's conditional
as a constant. On gaussian data the mode-based routes are measured
equivalent; on non-gaussian families the marginal is the robust one, the
measured comparison on a Poisson panel recovering the per-group
positions where the smoothed modes lose them.

## References

Muggeo, V. M. R. (2003). Estimating regression models with unknown
break-points. *Statistics in Medicine*, 22(19), 3055–3071.

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md),
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md),
[`seg_step()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))
dd$y <- 1 + 0.3 * dd$x + 2 * (dd$x > 6) +
  1.5 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)

built <- term_build(jseg(x), dd)
term_coef_names(built)
#> [1] "jseg.beta"   "jseg.gamma1" "jseg.delta1" "jseg.g1"    


# Fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(3)
  fd <- data.frame(x = sort(runif(300, 0, 10)))
  fd$y <- 1 + 0.3 * fd$x + 1.5 * (fd$x > 6) + 1.2 * pmax(fd$x - 6, 0) +
    rnorm(300, sd = 0.3)
  ft <- statmodels7::statmod(y ~ jseg(x),
                             distributions7::gaussian1_distrib(), fd)
  # truth: intercept 1, slope 0.3, and at 6 a jump of 1.5 with a change
  # of slope of 1.2
  round(coef(ft)$mu, 2)
}
#> (Intercept)   jseg.beta jseg.gamma1 jseg.delta1   jseg.psi1 
#>        1.08        0.28        1.19        1.62        6.00 
```
