# Segmented, Stepmented and Segmented-with-Jump Terms

A covariate whose effect changes at break-points estimated with
everything else. `seg` changes slope at each break-point and stays
continuous (muggeo2003); `jump` steps to a new level and is
discontinuous (fasola2018); `jseg` does both at the same points.

## Usage

``` r
seg(
  x,
  ...,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  penalty = NULL,
  c0 = 0.05,
  label = "seg"
)

jump(
  x,
  ...,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  penalty = NULL,
  c0 = 0.05,
  label = "jump"
)

jseg(
  x,
  ...,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  penalty = NULL,
  c0 = 0.05,
  label = "jseg"
)
```

## Arguments

- x:

  The covariate, an expression evaluated in the data.

- ...:

  At most one two-sided formula `psi ~ f` developing the break-points
  with covariates; see the section above. Cannot be combined with `by`.

- npsi:

  The number of break-points. Defaults to 1.

- psi:

  Optional starting positions; defaults to evenly spaced quantiles of
  the covariate. With a subformula they seed the development, each
  starting vector solving \\Z\gamma_k \approx \psi_k^{0}\\.

- by:

  An optional factor giving an independent set of break-points per
  level.

- linear:

  Whether the block carries the linear effect. Defaults to `TRUE`.

- penalty:

  The penalty on the changes: `NULL` (default, none), a penalties7
  penalty over as many coefficients as there are changes, or a function
  of that count returning one – a penalties7 constructor passed bare
  works, e.g. `penalty = penalties7::lasso_penalty`. A joint term
  declares two penalties, one on the slope changes and one on the jumps,
  and a penalty given as an object is used for both.

- c0:

  For a discontinuous term, the starting value of the scaling factor
  that separates the observations from the break-point, as a fraction of
  the distance to the ends of the range. Defaults to `0.05`, the value
  fasola2018 recommend; see Details.

- label:

  A single non-empty string prefixed to the coefficient names.

## Value

An object of class
[`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

Both constructions rest on the same device: a quantity that depends on
the break-point non-linearly is written as a linear function of it once
something is frozen at the previous iterate, so that a linear fit
returns the new break-point. Iterating that is the estimation algorithm,
and here it is the
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
contract of
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md) with a
different block: refreshing the term at the current coefficients and
fitting is one step of it.

### The continuous case

With \\f(x) = \delta (x-\psi)\_+\\ the contribution is differentiable in
\\\psi\\ away from the break-point, and \\\partial f/\partial\psi =
-\delta\\\mathbb{1}(x\>\psi)\\, so the design block is the ordinary
Jacobian and the break-point is an ordinary coefficient, updated by the
increment a Gauss-Newton step solves for. That is the algorithm of
muggeo2003 written in the coordinates the rest of this package uses: his
working variables \\U\\ and \\V\\ are its columns, and his update \\\psi
\leftarrow \psi + \gamma/\delta\\ is that step.

### The discontinuous case

A jump is not differentiable in \\\psi\\ at all, and yet needs no search
over candidate positions. The identity \$\$\mathbb{1}(x\>\psi) =
\frac{1}{2} + \frac{x-\psi}{2\lvert x-\psi \rvert}\$\$ holds exactly for
\\x \neq \psi\\, and is linear in \\\psi\\ once the weight \\1/(2\lvert
x - \psi\rvert)\\ is held at the previous iterate. Writing \\W =
1/(2\lvert x-\psi^{0}\rvert)\\ and \\Z = xW + 1/2\\, a jump of size
\\\kappa\\ at \\\psi\\ is \\\kappa Z + gW\\ with \\g = -\kappa\psi\\, so
a linear fit on \\(Z, W)\\ returns the break-point as \\\psi =
-g/\kappa\\ (fasola2018). The break-point is therefore not a coefficient
here but a quantity read off two of them, which is why refreshing the
term recovers it before rebuilding the weights.

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
covariates are computed on the rescaled covariate, while the truncated
line, the linear column and the reported contribution stay on the
original one.

The factor is not a constant. It governs how far the break-point may
travel in one step, so a large value lets the estimate leave a spurious
optimum and a small one is faithful to the step function. `c0` is its
starting value, and
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
halves it whenever the break-point reverses direction, which is the
signal that the iteration has begun to circle an optimum rather than
travel towards one. The run has converged when the change in every
break-point falls below a hundredth of the smallest distance between
distinct observations, which
[`seg_converged`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
reports.

### What the term carries

With `linear = TRUE` the block carries the linear effect too, so the
term is the whole relationship rather than the change in it. `by` gives
an independent set of break-points and changes per level of a factor.
`penalty` puts a penalty on the changes themselves – the slope changes
for `seg`, the jump sizes for `jump`, both for `jseg` – and leaves the
linear effect and the break-points alone; with the lasso that is a
selection of how many break-points are really there.

The penalty is declared through
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
which names the coefficients it covers, rather than attached to the
whole block through a map that selects them. The two describe the same
function of the same coefficients and are not interchangeable to a
fitting layer: a separable penalty under a selection map is the
generalized-lasso problem, whose proximal operator does not split by
coordinate, so
[`has_prox`](https://statmodels7.github.io/penalties7/reference/has_prox.html)
is `FALSE` for it and neither a proximal step nor a coordinate descent
can be taken. Named as coordinates the map is the identity and both are
available unchanged. `jseg` declares two penalties, one over the slope
changes and one over the jump sizes, since a change of slope and a
change of level are not comparable quantities and cannot share a
hyperparameter.

A break-point is confined to the interval between the 5th and the 95th
percentile of the covariate. Outside it the block is singular rather
than merely ill-conditioned: with \\\psi\\ below the smallest
observation the indicator is constant, so the truncated line and that
constant are linearly dependent, and with the linear effect present so
are all three columns. A run that ends against the limit has not located
a break-point, and
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
then returns the limit itself.

#### Break-points developed with covariates

A two-sided formula `psi ~ f` in `...` develops every break-point as
\\\psi_k = Z\gamma_k\\, with \\Z\\ the design of the right-hand side
through
[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
and one coefficient vector per break-point: `psi ~ g` with a factor `g`
is a break-point per group with the slopes shared, where `by` would give
every level its own slopes as well, and the two are therefore not
combinable. Each observation carries the position its own row of \\Z\\
implies, confined to the same interval as above.

For the continuous construction \\\gamma_k\\ are ordinary coefficients
and the block carries \\-\delta_k\\\mathbb{1}(x\>\psi_k)\\Z_j\\ in place
of the single Jacobian column, so a sub-term's penalty passes through
unchanged: `seg(x, psi ~ random(~1 | id))` is the random-changepoint
model of Muggeo, Atkins, Gallop and Dimidjian (2014), the per-group
positions shrunk towards a population one. For `jump` the identity above
splits the \\gW\\ column into \\W Z_j\\, whose coefficients are \\c_k =
-\kappa_k\gamma_k\\, and the development is read off as \\\gamma_k =
-c_k/\kappa_k\\ exactly as the scalar break-point is; a sub-term
carrying a penalty is rejected there, since the penalty would act on
\\c_k\\, which is the development scaled by the jump size, rather than
on \\\gamma_k\\ itself. `jseg` rejects a development: its reading of the
break-point is a quadratic in the increment that couples the slope
change with the jump, and it does not split over the columns of a
development, while the componentwise reading that remains diverges
whenever the jump size passes near zero mid-iteration.

The objective has local optima in the break-points, and the scaling
schedule widens the basin the iteration converges from rather than
removing the problem. Where the run begins therefore decides what it
finds, and
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
is the answer: it scores an equally spaced grid on the least-squares
profile and returns the specification with `psi` set to the best of it,
which is what fasola2018 recommend and what measurement supports over
both a conventional single start and bootstrap restarting. A continuous
term has no scaling factor, its working block being bounded already;
where its iteration alternates between two values the remedy is to
shrink the increment, as `segmented`'s `h` does.

## References

Muggeo, V. M. R. (2003). Estimating regression models with unknown
break-points. *Statistics in Medicine*, 22(19), 3055–3071.

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md),
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.md),
[`seg_step`](https://statmodels7.github.io/modelterms7/reference/seg_step.md),
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))
dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)

built <- term_build(seg(x), dd)
term_coef_names(built)
#> [1] "seg.lin"    "seg.delta1" "seg.psi1"  
seg_psi(built)
#> [1] 5.054907
```
