# Stepmented Term: a Level that Changes at Estimated Break-Points

A covariate whose effect is a constant that steps to a new value at
\\K\\ break-points estimated with everything else (fasola2018). Written
in the equation of any parameter of any distribution, the term
contributes

\$\$\sum\_{k=1}^{K} \delta_k\\\mathbb{1}(x_i \geq \psi_k)\$\$

to that parameter's linear predictor, so that in a model \\g(\theta_i) =
\eta_i\\ with the rest of the equation supplying \\z_i'\alpha\\,

\$\$\eta_i = z_i'\alpha + \sum\_{k=1}^{K} \delta_k\\\mathbb{1}(x_i \geq
\psi_k).\$\$

There is no linear effect: the relationship between \\x\\ and the
predictor is a step function, the intercept of the equation being its
level before the first break-point and \\\delta_k\\ the change of level
at \\\psi_k\\. A model that is linear in \\x\\ *and* steps is written by
putting the covariate in the equation as well, `y ~ x + jump(x)`, or by
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
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
  the covariate. Where a break-point carries a development they seed it,
  each starting vector solving \\Wp_k \approx \psi_k^{0}\\.

- by:

  A one-sided formula giving every coefficient of the term the same
  development, e.g. `by = ~0 + g` for an independent set per level of a
  factor. A bare variable is rejected; write the formula.

- c0:

  The starting value of the scaling factor that separates the
  observations from the break-point, as a fraction of the distance to
  the ends of the range. Defaults to `0.05`, the value fasola2018
  recommend; see Details.

- smoothed:

  `NULL` (the default: the construction exactly as documented above) or
  a penalties7
  [`abs_smoother`](https://statmodels7.github.io/penalties7/reference/abs_smoother.html),
  e.g.
  [`penalties7::smooth_probit()`](https://statmodels7.github.io/penalties7/reference/smooth_probit.html).
  The smoother replaces the step and the hinge by their smooth versions
  – \\(1 + s'(u))/2\\ and \\(u + s(u))/2\\ – so every break-point
  becomes an ordinary parameter of a \\C^\infty\\ model: there is no
  working parametrization, no auxiliary coefficient and no scaling
  schedule (`c0` is ignored, with a message), the block is the true
  Jacobian and the term is fitted by Gauss-Newton like
  [`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md). A
  development of a break-point – `psi ~ random(~1 | id)`, a penalized
  one included – is then legal for every kind, the read-off that
  constrained the discontinuous constructions having gone. The
  smoother's width is resolved at build from the covariate's spacing
  (the median gap between distinct values, within groups where a
  break-point development supplies a partition) unless the object
  carries one, and is reported: it is the width of the transition, the
  bent-cable reading, and the smoothing bias it buys is confined to a
  window of that width (probit, quintic) or decays as \\c/(4\|u\|)\\
  (hyperbolic). The objective is still multimodal in the positions –
  smoothing rounds the local optima, it does not remove them – so the
  profile start and the `n_boot` restarts stay necessary; a smoothed fit
  from a bad start has been measured converging to an absurd local
  optimum while reporting success.

- marginal:

  Whether the break-points are latent variables per group, integrated
  out of the likelihood exactly; see the section below. Requires the
  subformula `psi ~ random(~1 | g)`, whose grouping carries the latents.
  Defaults to `FALSE`, the construction documented above.

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
[`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)),
or of class
[`MarginalBreakTerm`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md)
with `marginal = TRUE`.

## Details

### How the break-points are estimated

A step is not differentiable in \\\psi\\ at all, and yet needs no search
over candidate positions. The identity \$\$\mathbb{1}(x\>\psi) =
\frac{1}{2} + \frac{x-\psi}{2\lvert x-\psi \rvert}\$\$ holds exactly for
\\x \neq \psi\\ and is linear in \\\psi\\ once the weight \\1/(2\lvert
x-\psi^{0}\rvert)\\ is held at the previous iterate. Writing \\W =
1/(2\lvert x-\psi^{0}\rvert)\\ and \\Z = xW + 1/2\\, a step of size
\\\delta\\ at \\\psi\\ is \\\delta Z + gW\\ with \\g = -\delta\psi\\, so
a linear fit on \\(Z, W)\\ returns the break-point as \\\psi =
-g/\delta\\ (fasola2018). The break-point is therefore not a coefficient
here but a quantity read off two of them, which is why refreshing the
term recovers it before rebuilding the weights, and why
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
exists.

### The scaling schedule

The weight is unbounded as \\x\\ approaches \\\psi\\, and since \\Z -
\psi W\\ is the step itself, \\Z\\ is \\\psi W\\ plus a quantity of
order one: an unbounded \\W\\ makes the two columns numerically
collinear and drowns the signal the fit reads. The remedy of fasola2018
is to move the observations rather than to cap the weight. With a
scaling factor \\c\\ the two intervals \\\[x\_{(1)}, \psi\]\\ and
\\(\psi, x\_{(n)}\]\\ are mapped onto \$\$\[x\_{(1)},\\ \psi - c(\psi -
x\_{(1)})\] \quad\text{and}\quad (\psi + c(x\_{(n)} - \psi),\\
x\_{(n)}\],\$\$ which leaves a gap of relative width \\c\\ around
\\\psi\\ and bounds \\W\\ without altering the model: the working
covariates are computed on the rescaled covariate, while the reported
contribution stays on the original one.

The factor is not a constant. It governs how far the break-point may
travel in one step, so a large value lets the estimate leave a spurious
optimum and a small one is faithful to the step function. `c0` is its
starting value and
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
halves it whenever the break-point reverses direction, the signal that
the iteration has begun to circle an optimum rather than travel towards
one. The run has converged when every break-point moves less than a
hundredth of the distance between consecutive distinct observations,
which
[`seg_converged`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
reports.

A break-point is confined to the interval between the 5th and the 95th
percentile of the covariate, outside which the indicator is constant and
the block singular, and the starting positions are chosen on a grid by
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.md):
measured over eight samples and four starting positions, the fraction of
runs recovering the break-point is 0 to 0.5 from a single conventional
start and 1 from the grid.

## The coefficients

`delta1` ... `deltaK`, the changes of level, followed by `g1` ... `gK`,
the auxiliary coefficients of the identity above, from which \\\psi_k =
-g_k/\delta_k\\ is read. The break-points are *not* coefficients; ask
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md).

## Developing a coefficient with covariates

As for
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md), a
two-sided formula in `...` whose left side names `delta`, `deltaK` or
`psi`, `psiK` develops that coefficient on covariates, and `by = ~f`
gives every coefficient the same development. The right side goes
through
[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
and takes any term of this package.

What the discontinuous construction can carry is narrower than what the
continuous one can, and the reason is the read-off rather than the
model. The auxiliary coefficient is \\g_k = -\delta_k\psi_k\\, a product
of the two quantities, so it stays linear in the unknowns only when one
of the two factors is a single number:

- **the break-point developed**, the change of level a single number:
  exact for any design, \\g_k = -\delta_k p_k\\ being linear in \\p_k\\.
  `jump(x, psi ~ id)` is a break-point per subject with a shared step
  size.

- **a development on group indicators**, where the design has one column
  per group and each observation belongs to one: the product collapses
  group by group and the read-off is exact within each.
  `jump(x, by = ~0 + g)` is an independent step and break-point per
  level.

- **anything else** is rejected rather than approximated, the product of
  two developments needing the outer product of their designs and no
  unconstrained fit returning it as one.

A sub-term carrying a penalty is rejected on a developed break-point,
since the estimated coefficients are \\-\delta_k p_k\\, the development
scaled by the step size, and a penalty would act on that rather than on
the development. A penalty on the changes of level themselves is
`jump(x, npsi = 4, delta ~ 0 + lasso(~1))`, the `0 +` removing the
subformula's own unpenalized intercept, which would otherwise be the
same column twice.

## The marginal construction

`jump(x, psi ~ random(~1 | g), marginal = TRUE)` treats each break-point
as a latent variable per group, \\\psi\_{ik} = m_k + u\_{ik}\\ with
independent \\u\_{ik} \sim N(0, \tau_k^2)\\, and integrates them out of
the likelihood exactly. The conditional is constant on the product
partition of the intervals between a group's ordered covariate values –
\\(n_i+1)^K\\ cells – but the sum is never taken over the cells: the
side process \\S_t = \\k : \psi_k \le x\_{(t)}\\\\ is monotone on the
subset lattice with independent coordinates, so it is a hidden Markov
chain on the \\2^K\\ side patterns whose transition factors over the
coordinates, and the exact forward recursion costs \\n K 2^K\\: a
coordinate that flips at a step contributes its interval's prior mass as
the transition weight, one that never flips its upper-tail mass at the
end. Up to eight break-points are covered; what prices more is not the
recursion but the \\2^K\\ components a fitting layer evaluates the
family at.

The prior is part of the likelihood: \\(m_k, \tau_k, \delta_k)\\ are
ordinary parameters estimated by maximum likelihood, with no penalty, no
marginal criterion and no smoothing constant, and the term is structural
([`MarginalBreakTerm`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md),
the likelihood shape of the contract) rather than a design block. With
one break-point the prior may be any continuous distributions7
distribution with its location fixed at zero, through
`random(distrib = )`: the interval masses are differences of its cdf and
their derivatives come from the cdf surface built for the censored
likelihoods, with that surface's own caveats (closed for the gaussian
and, in location and scale, for the t). `smoothed` and `c0` do not apply
and are ignored with a message; the posterior of each group's
break-points is read by
[`term_latent`](https://statmodels7.github.io/modelterms7/reference/term_latent.md).
On the step model the marginal is the route that resolves what a
smoothed mode cannot – the conditional is a step in the position, so a
Laplace approximation has no curvature to read – and it is the robust
route on non-gaussian families.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md),
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.md),
[`seg_step`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)

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
```
