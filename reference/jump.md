# Step Term: a Level that Changes at Estimated Break-Points

`jump()` specifies a covariate whose effect is a constant that steps to
a new value at \\K\\ break-points estimated with everything else
(Fasola, Muggeo and Kuchenhoff, 2018). Written in the equation of any
parameter of any distribution, the term contributes

\$\$\sum\_{k=1}^{K} \delta_k\\\mathbb{1}(x_i \geq \psi_k)\$\$

to that parameter's linear predictor, so that in a model \\g(\theta_i) =
\eta_i\\ in which the rest of the equation supplies \\z_i'\alpha\\,

\$\$\eta_i = z_i'\alpha + \sum\_{k=1}^{K} \delta_k\\\mathbb{1}(x_i \geq
\psi_k).\$\$

There is no linear effect: the relationship between \\x\\ and the
predictor is a step function; the intercept of the equation is its level
before the first break-point, and \\\delta_k\\ is the change of level at
\\\psi_k\\. A model that is linear in \\x\\ *and* steps is written by
putting the covariate in the equation as well, `y ~ x + jump(x)`, or by
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
where the slope is to change at the same points.

## Usage

``` r
jump(
  x,
  ...,
  npsi = 1,
  psi = NULL,
  by = NULL,
  c0 = 0.05,
  smoothed = NULL,
  marginal = FALSE,
  n_boot = 10,
  label = "jump"
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
  the covariate. Where a break-point carries a development, they seed
  it, and each starting vector solves \\Wp_k \approx \psi_k^{0}\\.

- by:

  A one-sided formula giving every coefficient of the term the same
  development, e.g. `by = ~0 + g` for an independent set per level of a
  factor. A bare variable is rejected; write the formula.

- c0:

  The starting value of the scaling factor that separates the
  observations from the break-point, as a fraction of the distance to
  the ends of the range. Defaults to `0.05`, the value that Fasola et
  al. (2018) recommend; see Details.

- smoothed:

  `NULL` (the default: the construction exactly as documented above) or
  a numericals7
  [`numericals7::abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.html),
  e.g.
  [`numericals7::smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.html).
  The smoother replaces the step and the hinge by their smooth versions,
  \\(1 + s'(u))/2\\ and \\(u + s(u))/2\\, so every break-point becomes
  an ordinary parameter of a model as smooth as the smoother
  (\\C^\infty\\ for the probit and the hyperbolic, \\C^3\\ for the
  quintic; the page of the quintic gives its cost under an outer
  criterion, where the probit is the recommended smoother). There is no
  working parametrization, no auxiliary coefficient and no scaling
  schedule (for `jump()` and
  [`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
  a `c0` that is passed is ignored with a message). The block is the
  true Jacobian and the term is fitted by Gauss-Newton like
  [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md). A
  smoothed term accepts a development of any coefficient, including a
  break-point (`psi ~ random(~1 | id)`) and penalized developments, for
  every kind, since a smoothed term has no read-off that would constrain
  the development. The smoother's width is resolved at build from the
  covariate's spacing (the median gap between distinct values, within
  groups where a break-point development supplies a partition) unless
  the object carries one. For the quintic, which is exact outside its
  width, the width is raised to at least 0.55 times the largest gap over
  the range that a break-point may take, so that a break-point cannot
  sit in a gap with no observation within its width. The width is
  reported: it is the width of the transition (the bent-cable reading).
  The smoothing bias is zero outside the transition for the quintic,
  decays with the Gaussian tails for the probit, and decays as
  \\c/(4\|u\|)\\ for the hyperbolic, whose parameter \\c\\ is a squared
  length. The objective remains multimodal in the positions: smoothing
  rounds the local optima and does not remove them, so the profile start
  and the `n_boot` restarts remain necessary.

- marginal:

  Whether the break-points are latent variables per group, integrated
  out of the likelihood exactly; see the section below. Requires the
  subformula `psi ~ random(~1 | g)`, whose grouping carries the latents.
  Defaults to `FALSE`, the construction documented above.

- n_boot:

  How many bootstrap restarts the fitting layer runs after the iteration
  first converges (Wood 2001, the device that `segmented` runs by
  default): each restart re-estimates on a bootstrap resample from the
  current break-points and then on the data again from where the
  resample ended, keeping the better fit. The objective has local optima
  in the break-points, and with several of them a grid start does not
  reach the right basin from every placement. Defaults to 10,
  `segmented`'s own default; 0 disables. The term itself only declares
  the value; running the restarts is the fitting layer's, as with a
  penalty's hyperparameters.

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

### How the break-points are estimated

A step is not differentiable in \\\psi\\, and the identity
\$\$\mathbb{1}(x\>\psi) = \frac{1}{2} + \frac{x-\psi}{2\lvert x-\psi
\rvert}\$\$ holds exactly for \\x \neq \psi\\ and is linear in \\\psi\\
once the weight \\1/(2\lvert x-\psi^{0}\rvert)\\ is held at the previous
iterate, so a search over candidate positions is unnecessary. Writing
\\W = 1/(2\lvert x-\psi^{0}\rvert)\\ and \\Z = xW + 1/2\\, a step of
size \\\delta\\ at \\\psi\\ is \\\delta Z + gW\\ with \\g =
-\delta\psi\\, so a linear fit on \\(Z, W)\\ returns the break-point as
\\\psi = -g/\delta\\ (Fasola et al., 2018). The break-point is therefore
not a coefficient here but a quantity read off two of them, which is why
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
recovers it before rebuilding the weights, and why
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
exists.

### The scaling schedule

The weight is unbounded as \\x\\ approaches \\\psi\\, and since \\Z -
\psi W\\ is the step itself, \\Z\\ is \\\psi W\\ plus a quantity of
order one: an unbounded \\W\\ makes the two columns numerically
collinear and drowns the signal that the fit reads. The remedy of Fasola
et al. (2018) is to move the observations and not to cap the weight.
With a scaling factor \\c\\ the two intervals \\\[x\_{(1)}, \psi\]\\ and
\\(\psi, x\_{(n)}\]\\ are mapped onto \$\$\[x\_{(1)},\\ \psi - c(\psi -
x\_{(1)})\] \quad\text{and}\quad (\psi + c(x\_{(n)} - \psi),\\
x\_{(n)}\],\$\$ which leaves a gap of relative width \\c\\ around
\\\psi\\ and bounds \\W\\ without altering the model: the working
covariates are computed on the rescaled covariate, while the reported
contribution stays on the original one.

The factor changes during the iteration and governs how far the
break-point may travel in one step, so a large value lets the estimate
leave a spurious optimum and a small one is faithful to the step
function. `c0` is its starting value and
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
halves it whenever the break-point reverses direction, the signal that
the iteration has begun to circle an optimum instead of traveling toward
one. The run has converged when every break-point moves less than a
hundredth of the median distance between consecutive distinct
observations, which
[`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
reports (see
[`seg_step()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
for the relation to the smallest distance that Fasola et al. use).

A break-point is confined to the interval between the 5th and the 95th
percentile of the covariate, so that both sides of a break-point carry
observations (a break-point beyond the range of the data would make the
indicator constant and the block singular). The starting positions are
chosen on a grid by
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md).

## The coefficients

The coefficients are `delta1` ... `deltaK`, the changes of level,
followed by `g1` ... `gK`, the auxiliary coefficients of the identity
above, from which \\\psi_k = -g_k/\delta_k\\ is read. For a sharp term
the break-points are not coefficients, and
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
reports them. A smoothed term holds them directly, and its coefficients
are `delta1` ... `deltaK` followed by `psi1` ... `psiK`.

## Developing a coefficient with covariates

As for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md), a
two-sided formula in `...` whose left side names `delta`, `deltaK` or
`psi`, `psiK` develops that coefficient on covariates, and `by = ~f`
gives every coefficient the same development. The right side goes
through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
and takes any term with a fixed design block; structural terms and terms
whose block moves with its coefficients are rejected.

The sharp discontinuous construction can carry a narrower set of
developments than the continuous one, because of the read-off. The
auxiliary coefficient is \\g_k = -\delta_k\psi_k\\, a product of the two
quantities, so it stays linear in the unknowns only when one of the two
factors is a single number:

- The break-point is developed and the change of level is a single
  number: the read-off is exact for any design, because \\g_k =
  -\delta_k p_k\\ is linear in \\p_k\\. `jump(x, psi ~ id)` is a
  break-point per subject with a shared step size.

- A development on group indicators, where the design has one column per
  group and each observation belongs to one: the product collapses group
  by group and the read-off is exact within each. `jump(x, by = ~0 + g)`
  is an independent step and break-point per level.

- Anything else is rejected and is not approximated, since the product
  of two developments needs the outer product of their designs, and an
  unconstrained fit cannot return that as one set of coefficients.

A smoothed term accepts a development of any coefficient.

A sub-term carrying a penalty is rejected on a developed break-point of
a sharp term, since the estimated coefficients are \\-\delta_k p_k\\,
the development scaled by the step size, and a penalty would act on that
product and not on the development. A penalty on the changes of level
themselves is `jump(x, npsi = 4, delta ~ 0 + lasso(~1))`. The `0 +`
removes the subformula's own unpenalized intercept, which would
otherwise repeat the constant column of the sub-term, a design that is
rejected when the term is built.

## The marginal construction

`jump(x, psi ~ random(~1 | g), marginal = TRUE)` treats each break-point
as a latent variable per group, \\\psi\_{ik} = m_k + u\_{ik}\\ with
independent \\u\_{ik} \sim N(0, \tau_k^2)\\, and integrates them out of
the likelihood exactly. The parameters \\m_k\\ and \\\tau_k\\ are named
`psik.mean` and `psik.sd`. The conditional is constant on the product
partition of the intervals between a group's ordered covariate values,
giving \\(n_i+1)^K\\ cells. The sum is never taken over the cells. The
side process \\S_t = \\k : \psi_k \le x\_{(t)}\\\\ is monotone on the
subset lattice with independent coordinates, so it is a hidden Markov
chain on the \\2^K\\ side patterns whose transition factors over the
coordinates, and the exact forward recursion costs \\n K 2^K\\. In this
recursion a coordinate that flips at a step contributes its interval's
prior mass as the transition weight, and a coordinate that never flips
contributes its upper-tail mass at the end. Up to eight break-points are
covered. The cost of more lies in the \\2^K\\ components at which a
fitting layer evaluates the family.

The prior is part of the likelihood: \\(m_k, \tau_k, \delta_k)\\, named
`psik.mean`, `psik.sd` and `deltak`, are ordinary parameters estimated
by maximum likelihood, with no penalty, no marginal criterion and no
smoothing constant, and the term is itself structural
([`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md),
of the likelihood shape). With one break-point the prior may be any
continuous distributions7 distribution whose location is fixed at zero,
through `random(distrib = )`: the interval masses are differences of its
cdf and their derivatives come from the cdf surface built for the
censored likelihoods, with that surface's own caveats (closed for the
gaussian and, in location and scale, for the t). `smoothed`, `c0` and
`n_boot` do not apply and are ignored with a message; the posterior of
each group's break-points is read by
[`term_latent()`](https://statmodels7.github.io/modelterms7/reference/term_latent.md).
The marginal construction does not use a Laplace approximation of the
latent break-points, which the mode-based routes need. On the step model
the conditional is a step in the position, so such an approximation has
no curvature to read.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md),
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md),
[`seg_step()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))
dd$y <- 1 + 2 * (dd$x > 6) + rnorm(200, sd = 0.3)

built <- term_build(jump(x), dd)
term_coef_names(built)
#> [1] "jump.delta1" "jump.g1"    
seg_psi(built)
#> [1] 5.054907


# The term is fitted next. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(2)
  fd <- data.frame(x = sort(runif(200, 0, 10)))
  fd$y <- 1 + 2 * (fd$x > 6) + rnorm(200, sd = 0.3)
  ft <- statmodels7::statmod(y ~ jump(x),
                             distributions7::gaussian1_distrib(), fd)
  # The truth is intercept 1, a jump of 2 and the break-point at 6.
  round(coef(ft)$mu, 2)
}
#> (Intercept) jump.delta1   jump.psi1 
#>        1.02        2.01        5.96 
```
