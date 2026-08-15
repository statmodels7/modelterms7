# modelterms7 0.44.0

* `s()` and `te()` take `sparse`, and a FACTOR `by` is where it applies: each
  row sits in the block of its own level and nowhere else, a density of `1/m`,
  which is the shape `random()`'s block already had. Measured at 2000 rows,
  `k = 10` and 200 levels, 0.35 MB against 28.93 MB, the numbers and the
  coefficient names identical, and the storage is in the blueprint so
  prediction does not densify.

* It is REFUSED where there is nothing to build on rather than accepted and
  ignored: without a `by` the basis is dense by construction, the
  Demmler-Reinsch rotation making it so, and a numeric `by` merely multiplies
  it. The first is caught at construction and the second at the build, which
  is where each becomes knowable.

* The BLOCK alone is sparse. The penalty of a factor `by` is the same matrix
  repeated blockwise and comes back dense -- 25.92 MB at those sizes, at a
  density of 0.0005 -- because `penalties7` densifies at its contract
  boundary. Making that sparse is a change to that package's contract, not to
  this construction, and it is not made here.

# modelterms7 0.43.0

* `linpar()` and the FORMULA route of the five penalized terms take
  `sparse = FALSE`. It builds the block through
  `Matrix::sparse.model.matrix()`, which BUILDS it sparse rather than building
  a dense matrix and compressing it: measured at 20000 rows and a factor of
  1000 levels, 0.002 s and 1.8 MB against `stats::model.matrix`'s 0.100 s and
  161.5 MB, the numbers identical, and a design that would be 32 GB dense
  builds in 0.02 s and 19 MB -- which is what says there is no dense
  intermediate.

* It pays where the formula carries a factor of many levels, whose indicator
  columns hold one non-zero per row, and costs more than it saves on numeric
  covariates, whose block is dense whatever is asked for. The storage is part
  of the blueprint, so `term_predict()` builds new data the same way, and
  standardization does not undo it, being a diagonal map on the penalty and
  never an operation on the design.

* `random()` does NOT take it: `.random_block()` builds with
  `Matrix::sparseMatrix()` unconditionally and its density is `1/m` by
  construction, so the argument could only have densified. A matrix input to
  any of the five does not take it either, being kept in whatever storage it
  arrives in.

* `linpar()` takes `contrasts`, the coding for its factors, which was
  reachable only through the session's `options("contrasts")`. The blueprint
  carries it, so prediction is coded the way the fit was.

* `interpret_formula()` takes `linpar`, arguments for the IMPLICIT linpar --
  the term the bare covariates collapse into, which a caller never writes, so
  it is the only place they can be given.

# modelterms7 0.42.0

* `enet()`, `scad()` and `mcp()` take `search`, how their own two
  hyperparameters are combined: `"grid"` for every combination of them,
  `"cyclic"` for one at a time. It sat on the criterion, where it did not
  belong -- the same criterion is put to the smooth hyperparameters of a
  model, which are read at the mode rather than swept, so most of what it is
  asked about could not use it. A penalty with a kink is fitted by a scheme
  of its own, and how that scheme covers the term's own hyperparameters is
  part of the scheme. `term_search()` reports it.

* `lasso()` and `ridge()` do not carry it: with one hyperparameter or none
  there is nothing to combine, and the argument is reported by name rather
  than accepted and read by nothing.

* `by` is gone from all five, and it was not merely unimplemented: there is
  nothing for it to build. Where `by` earns its place, in `s()` and `te()`,
  it builds PENALTY STRUCTURE -- one block of the penalty matrix per level of
  the factor. These five are separable, which is the property that lets a
  coordinate descent fit them, so repeating the penalty blockwise is the
  identity operation and `lasso(~ (x1 + x2):g)` already IS the per-level fit,
  coefficient by coefficient. It was reserved for a later release that had no
  work to do.

* `ridge()` loses `n_lambda` and `min_ratio`, which it never read. They were
  accepted "so that the five constructors read alike", and a symmetry of
  spelling is not worth an argument that does nothing: a ridge is twice
  differentiable everywhere, so its `lambda` is estimated at the mode by
  `reml()` and there is no path for a grid size or a depth to describe. Both
  are now reported by name, as any other argument the term does not have.

* The defaults are NUMBERS ON THE SIGNATURE: `n_lambda = 25`, `n_alpha = 5`,
  `n_a = 5`, `n_gamma = 5`, `min_ratio = 1e-4`, `search = "grid"`. They were
  `NULL` and documented as "leaves it to the criterion", which named a number
  nothing printed -- a default a reader cannot find is a default a reader
  cannot choose against. The two grid sizes differ because the axes do: a
  path over the size of the kink spans four decades and wants that many
  points, an axis beside it spans one bounded interval and does not.

# modelterms7 0.40.0

* A multivariate Student t prior now carries standard errors and intervals
  for every hyperparameter, its family having supplied the third and fourth
  derivatives in the response. The page no longer names it as an exception:
  where an interval is absent the cause is the POINT and not a missing
  derivative.

# modelterms7 0.39.0

* A location HELD away from zero is no longer rejected. What is unidentified
  is a FREE location, which competes with the intercept of the equation the
  term sits in; a held one shrinks the effects towards its value, which is a
  modelling statement. Nor could the value be policed in general: where the
  prior is a transformation of another family the parameter is the mean on
  the ORIGINAL scale, and holding a gamma's mean at one is what centers its
  logarithm -- at a value zero would put outside the parameter's own domain.
  `fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)` is now a
  usable effects distribution, and its free-location form is still refused.

* The page says which hyperparameters carry an interval and which do not, and
  why: a run that ended where the criterion has no maximum has a curvature of
  the wrong sign, and a multivariate Student t prior has no exact outer
  curvature until that family supplies its third and fourth derivatives in
  the response.

# modelterms7 0.38.0

* A multivariate Student t is admitted as the effects' distribution, its
  mixed response-parameter block having been written. The admissibility rule
  is unchanged and was already a property rather than a list of names, so
  nothing here had to learn about the family: it carries a location block as
  long as its dimension and a scale matrix, and it answers the block a
  marginal criterion reads.

# modelterms7 0.37.0

* `random()` says two things and nothing else: the formula, and the
  distribution of the effects. Which chart the hyperparameters ride, what
  they are called, how many there are and where the log-density has a kink
  are properties of that distribution, read off it rather than restated.

* The hyperparameter of `random(~ 1 | g)` IS the standard deviation of the
  effects. It was `lambda`, a precision, which a reader of a mixed model had
  to invert and take the root of; a variance component is now reported as
  one. It remains the same model as a ridge, which a test pins.

* `distrib` accepts a MULTIVARIATE distribution of the within-group
  dimension, which is what lets the effects of a group depend on each other:
  `mvgaussian_distrib(2, omega = ar1(2))` is a prior whose precision is
  autoregressive, and the free names then say WHICH matrix the structure is.
  Correlation is admitted by a PROPERTY -- a location block as long as the
  dimension, and a covariance, precision or scale matrix -- so a family added
  later is covered and the simplex-valued ones are refused without being
  named. A family with no mixed response-parameter block is rejected at the
  build rather than at the criterion, where the message would name a generic.

* A UNIVARIATE distribution over several within-group columns is a TEMPLATE:
  one copy per column with its own hyperparameters, declared as one
  `term_penalties()` entry per column over that column's own stride. An
  intercept and a slope are quantities of different units, and one shared
  scale would price them against each other -- with a pseudo-Huber it would
  share the point where the loss stops being quadratic as well. A list of
  distributions gives one per column explicitly.

* A prior on the effects is CENTERED: a free location is confounded with the
  intercept of the equation the term sits in, so it is rejected with the
  spelling that fixes it rather than fitted along a flat direction.

* `precision` is removed. A structured precision is the matrix parameter of a
  multivariate gaussian, so it is written as one -- which also says which of
  the two matrices the structure is. The argument is reported by name rather
  than swallowed by the dots.

* `kinks` is removed, and it was not merely redundant. Its default of
  `numeric(0)` OVERRODE `penalties7::distrib_kinks()`, so a Laplace prior on
  a random effect declared no kink where a `lasso()` term on the same parent
  declares one at zero -- and a fitting layer chooses the scheme by asking,
  so the block went to the system solved by a curvature the penalty does not
  have there. `NULL`, the spelling that would have fixed it, was rejected by
  the property's own type.

* `correlated` beside `distrib` is an error rather than ignored; the page
  said it was ignored, which is worse than an error.

# modelterms7 0.36.0

* A hyperparameter's argument carries a third state: several numbers are the
  grid a path visits, used as they stand. `lasso(~x, lambda = c(0.1, 1, 10))`
  sweeps exactly those three, and `enet(~x, lambda = seq(0.1, 10, length = 10))`
  writes out the grid of `lambda` while leaving `alpha` to be built, the state
  being settled per hyperparameter and not per term. `term_values()` reports
  them, beside `term_hyper()` for the ones held at one number.

* A written-out grid is still ESTIMATED: what the caller fixed is where to
  look, not the answer, so it is not reported as held. The value that empties
  the block does not cap it and `min_ratio` does not extend it -- both
  construct a grid, and here there is nothing to construct.

* Several values are rejected wherever the penalty has no kink, rather than
  taken and ignored: there is no path to visit them on, the hyperparameter
  being read at the mode by the criterion. The question is put to the penalty
  at a probe value and not written per constructor, so `ridge()` and a
  `random()` effect under a Gaussian prior are covered by one line while the
  same effect under a Laplace prior is not.

# modelterms7 0.35.0

* `ridge()`, `lasso()`, `enet()`, `scad()` and `mcp()` have a page each,
  carrying that penalty's own formula, its hyperparameters and the interval
  each of them may lie in. What the five share -- the formula and matrix
  input, the standardization, the prediction -- is documented once on
  `penalized_terms`, which they inherit their common arguments from.

* Each takes the SIZE OF THE GRID for each of its hyperparameters:
  `lasso(~x, n_lambda = 50)`, `enet(~x, n_lambda = 40, n_alpha = 12)`,
  `scad(~x, n_a = 6)`, `mcp(~x, n_gamma = 7)`. How finely a hyperparameter
  is swept belongs to the term for the same reason as whether it is swept
  at all: a penalized block of four columns and one of four hundred want
  different grids, and a criterion applies to every term of the model at
  once and cannot know which it is looking at. `NULL`, the default, leaves
  it to the criterion. `term_grid()` reports it, and the value travels with
  the penalty's entry as the held values do.

* And `min_ratio`, how far down the path reaches as a fraction of the kink
  that empties the block. One number per term rather than one per
  hyperparameter, because only the sweep by kink size uses it: a bounded
  hyperparameter is swept over its own interval and a shape that does not
  move the kink over a geometric grid above its lower bound, where a
  fraction of an emptying value means nothing. `term_path_min()` reports it.

# modelterms7 0.33.0

* `ridge()` takes `lambda`, the PRECISION of the Gaussian prior, where it
  took `sigma`: larger means more shrinkage, as it already did for
  `lasso()`, and it is the same number the quadratic penalty behind `s()`
  and `te()` calls by that name. `random()`'s default penalty follows, so
  the variance component a reader wants from it is `1/sqrt(lambda)`.

# modelterms7 0.32.0

* Every term takes its penalty's hyperparameters as arguments, and holds
  the ones it is given: `lasso(x, lambda = 3)`, `ridge(x, sigma = 0.5)`,
  `enet(x, lambda = 2, alpha = 0.5)`, `scad(x, lambda = 1, a = 3.7)`,
  `mcp(x, gamma = 3)`, `s(x, lambda = 2)`, `te(x, z, lambda = c(1, 5))`,
  `random(~1 | g, hyper = c(sigma = 0.4))`. NULL, the default, means the
  hyperparameter is estimated. Which ones are held is a property of the
  term, since the term is where the penalty is named.

  The arguments carry the penalty's OWN names -- `sigma` for a Gaussian
  prior, `a` for SCAD -- which are the names a summary prints and
  penalties7 documents. `ridge(x, lambda = 2)` is rejected with a message
  naming `sigma`, where R would have reported an unused argument.

* `term_hyper()` reports what a term holds, and every entry of
  `term_penalties()` carries its own held values in the field `fixed`, so a
  structural term propagates what its sub-terms hold by copying the entry:
  `gas(by = ~ ridge(id, sigma = 2))` and `nl(a ~ 0 + ridge(~g, sigma = 0.5))`
  reach the fit with the value the caller wrote.

# modelterms7 0.31.0

* `term_coef_start()`, a new generic: the coefficients a built term asks to
  be started at. The base method returns zero everywhere, which is what a
  term whose block is a fixed design wants, and a term that recomputes its
  block from its coefficients answers with the start `term_build()`
  computed. Zero is degenerate rather than neutral there: in `jump()` the
  break-point is read off two coefficients as `-g_k/delta_k`, so a vector
  of zeros puts every break-point at the same clamped position and leaves
  the block singular, and in `seg()` the Jacobian column vanishes
  identically. `nl()` answers with the values its own `start` names.

* `term_readable()` for a break-point term: the quantities of the model it
  defines -- the linear effect, the changes of slope and of level, and the
  break-points `psi_k` -- with the Jacobian from the coefficients, so a
  caller holding their variance matrix carries it across. A continuous
  term holds each of them as a coefficient; a discontinuous one does not
  hold the break-point at all, and the rows

      d psi_k / d g_k = -1 / delta_k,   d psi_k / d delta_k = g_k / delta_k^2

  are the delta method `segmented` reports a break-point's standard error
  by. Checked against numDeriv on all three constructions (1e-11). A
  developed coefficient has no single number to report and the method
  answers `NULL`.

# modelterms7 0.30.0

* `seg()`, `jump()` and `jseg()` have a documentation page each, opening
  with the model the term defines. Written for the predictor of any
  parameter of any distribution rather than for a gaussian mean, they are

      seg    eta_i = z_i'alpha + beta x_i + sum_k gamma_k (x_i - psi_k)_+
      jump   eta_i = z_i'alpha            + sum_k delta_k 1(x_i >= psi_k)
      jseg   eta_i = z_i'alpha + beta x_i
                     + sum_k [delta_k 1(x_i >= psi_k) + gamma_k (x_i - psi_k)_+]

  with the rest of the equation supplying `z_i'alpha`.

* The coefficients are named after those formulas: `beta` for the linear
  effect (was `lin`), `gamma_k` for a change of slope (was `delta_k`) and
  `delta_k` for a change of level (was `kappa_k`). `psi_k` is unchanged for
  the continuous construction, and the auxiliary pair of the discontinuous
  ones stays `g_k`.

* `jump()` has no `linear` argument. A step model's relationship with the
  covariate is a step function, the intercept of the equation being its
  level before the first break-point; a model that is linear in the
  covariate and steps is `y ~ x + jump(x)`, and one whose slope changes at
  the same points is `jseg()`.

* `by` is a one-sided formula and develops EVERY coefficient of the term on
  it: `by = ~0 + g` is the independent set per level that `by = g` used to
  give. A bare variable is rejected, with the formula it stands for in the
  message.

* Any coefficient of the term takes a development, through a two-sided
  formula in `...` whose left side names it -- `psi`, `psi1`, `gamma`,
  `gamma2`, `delta`, `beta` -- a stem meaning every coefficient of its
  kind. The right side goes through `interpret_formula()` and takes any
  term of the package. The continuous construction accepts any combination;
  the discontinuous ones read the break-point off a product of the unknowns
  and accept a developed break-point alone (exact on any design), or every
  varying coefficient on one shared design with a column per group.

* `penalty=` is retired from all three. A penalty on the changes is the
  development on a penalized intercept, `seg(x, npsi = 4,
  gamma ~ 0 + lasso(~1))`, which is the lasso that selects how many
  break-points there are; the entries of a subformula shared by every
  coefficient of a kind are pooled into one penalized block under one
  hyperparameter.

* A penalized constructor keeps the intercept where it is all the formula
  has, so `lasso(~1)` is a block of that one column rather than a block of
  none. That is what makes the spelling above mean anything.

* `interpret_formula()` removes a bare covariate that a `seg()` or `jseg()`
  term in the same equation already carries the linear effect of, and
  reports the removal: the two are exactly the same column, so
  `y ~ x + seg(x)` was rank deficient by one. `seg(x, linear = FALSE)`
  keeps the linear effect outside the term instead. Another term spanning
  the same direction, as a spline basis does, is reported without being
  modified.

# modelterms7 0.29.0

* `term_third()`, a new generic: the second derivative of a structural
  term's predictor differentiated once more and contracted against ONE
  direction, together with the derivative of the term's jacobian along the
  same direction. It is what the exact gradient of a marginal criterion
  needs where a penalty covers the term's own parameters, the criterion
  asking for `tr(M dK/du[v])` at the direction the penalized mode moves in.
  The full third derivative is an `m^3` array per observation and is never
  formed: what is propagated is a matrix per observation, the same size as
  the curvature and therefore the same O(n m^2). Measured against
  `term_curvature()` on panels of 9, 15 and 30 unknowns: 2.62x, 2.66x and
  2.78x, flat in the number of unknowns.

  The base method on `model_term` returns zero, which is right for a term
  whose predictor is a block of columns; `structural_term` REFUSES, so a
  term that bends the predictor and has not written its third derivative
  reports nothing rather than a zero a caller could not tell from a genuine
  one. `GasTerm` implements it on both routes, scalar and submodel.

  Each order of differentiating the predictor through the recursion pulls in
  one more order of the family, the score the recursion is driven by being
  read at the predictor it produces: the curvature's `M` is built from third
  derivatives and this needs a FOURTH, which `blocks` supplies as `N`
  alongside `dcurv`. Validated against a central difference of
  `term_curvature()` along the direction, Richardson-extrapolated because
  the plain difference's own truncation is larger than the gap being
  measured: 4.5e-11 to 1.7e-9 relative over p and q to 3, on a stub whose
  four derivative orders are all non-zero and all bounded.

* `gas_levinson3()`: the Levinson-Durbin map's third derivative, contracted
  against one direction. ⚠️ It is identically zero for q <= 2, the map being
  multilinear of degree k in the first k partial autocorrelations, so a
  check that stops at q = 2 compares zero with zero; the test runs to q = 4.

* `.gas_curvature_core()`: the second and third orders share one body rather
  than a method each. The third order reads F, Phi and the score's first two
  derivatives at every lag, so a separate implementation would carry a
  second copy of the first two orders and the two would drift.

# modelterms7 0.28.0

* The general recursion of the submodel route runs compiled
  (`gas_filter_sub_cpp` in `src/gas_filter.cpp`). The R side prepares each
  group's derivative rows once, scattered onto the group's active
  coordinates by one vectorized assignment per parameter, so the kernel's
  loop allocates nothing per step; the score and curvature callbacks stay
  in R, as in the scalar kernel and for the same reason. Measured on a
  panel with a developed level: 4.2x-4.3x over the R route, flat in n
  (n = 2000 to 8000), the compiled route now within 1.8x-2.1x of the
  scalar kernel where the R route sat at 5.4x-5.6x. The R loop survives
  as `.gas_filter_sub_r`, the twin the kernel is compared against at a
  tolerance (a compiler may contract a multiply-add; measured agreement
  is one rounding, 2.2e-16 on the jacobian). The adjoint and the
  curvature stay in R: the adjoint is one pass computed once per
  gradient, and the curvature is computed once, at `vcov()`.

# modelterms7 0.27.0

* Every parameter of `gas()` takes a subformula: a two-sided formula in
  `...` develops it as `psi_{j,t} = g_j^-1(z_t' gamma_j)` over the design
  of the right-hand side through `interpret_formula()`, so
  `gas(p = 1, q = 1, omega ~ ridge(~g), alpha1 ~ s(x), pacf1 ~ random(~1 | id))`
  is a specification. The development acts on the unconstrained scale of
  the parameter's own link, so every per-observation value stays in its
  own set whatever the coefficients are -- a loading on the log link is
  positive at every observation, a persistence on the rhobit chart is
  stationary at every observation. The coefficients are the term's
  parameters, unconstrained on the identity link, with the chart applied
  inside; the sub-terms' penalties are reported through
  `term_penalties()` under `parameter::subterm`.

  A parameter that varies by observation changes the recursion itself,
  `f_t = omega_t + sum_i a_{i,t} s_{t-i} + sum_j b_{j,t} f_{t-j}`, with
  `b_t` from the Levinson-Durbin map of that observation's partial
  autocorrelations, so the filter, the reverse pass and the second-order
  propagation all run a general per-observation recursion (in R;
  `src/gas_filter.cpp` keeps the scalar case). Validated each piece
  against the special case it generalizes and against `numDeriv`: a
  development of `~1` reproduces the scalar filter to 1e-12 with the
  chart's chain factor exact; before the shorthand was retired,
  `omega ~ random(~1 | g)` with `by = g` was pinned against
  `deviations = TRUE` to 1e-12, column for column, on the identity chart
  and on the log one, and at the FIT level the two agreed on the REML
  hyperparameter to the printed digit; a time-varying loading and a
  per-group development each reproduce a recursion written by hand; and
  the jacobian, the adjoint and the curvature agree with `numDeriv` to
  the reference's own accuracy at p and q up to 2, a developed partial
  autocorrelation included.

  `by = ~f` (a formula, where a grouping variable is a bare symbol)
  gives every parameter the same subformula; mixing it with
  per-parameter formulas is an error. Everything is carried on each
  group's ACTIVE SET -- a development's coordinate reaches only the
  groups where its column is not identically zero -- so with grouping
  indicators the per-observation work is constant in the number of
  groups, the property the deviations machinery had and the measured
  twelve-minute trap of a full square at five hundred groups required.
  Measured at n up to 10000, the general recursion costs 5.4-5.6x the
  compiled scalar filter, flat in n; porting it is the natural next step
  and the scalar kernel's own measurement (callbacks 17-27 per cent of
  the loop) says roughly what it would pay.

* `gas(deviations =, penalty =)` and `nl(penalty =, penalize =)` are
  RETIRED: the subformulas subsume them (Giovanni, explicit). A gas
  panel's population-and-departures model is
  `gas(omega ~ random(~1 | id), by = id)`, proved equivalent to the
  shorthand to the printed digit before the removal; a penalized nl
  parameter is `a ~ lasso(~z1 + z2)` or `a ~ ridge(~g)`, the sub-term
  bringing its own hyperparameter. What the removal loses is only the
  penalty on a SCALAR nl parameter's single coefficient, which had no
  subformula spelling. `seg()`, `jump()` and `jseg()` KEEP `penalty=`:
  it reaches the changes (the slope changes and the jump sizes), which
  are not parameters a subformula models, and a lasso there is the
  selection of how many break-points are real. A call carrying a removed
  argument is reported by name.

* `penalty` takes only a penalties7 object or a function of the
  coefficient count -- the strings `"none"`/`"lasso"`/`"ridge"` are gone
  (Giovanni, explicit), `NULL` being the default. A function with an
  `n_coef` formal is called by that NAME, so a penalties7 constructor
  passes bare: `penalty = penalties7::lasso_penalty` is exactly as short
  as the string it replaces, where before it failed with "the function
  covers 1 coefficients" because the constructor's first formal is the
  map and received the count positionally.

* The break-points of `seg()` and `jump()` take a subformula,
  `psi ~ f`, developing every break-point as `psi_k = Z gamma_k` over the
  design of the right-hand side through `interpret_formula()`.

  `psi ~ g` with a factor is a break-point per group with the slopes
  shared (where `by` would give every level its own slopes as well, so
  the two are not combinable), and each observation carries the position
  its own row implies, confined to the same [q05, q95] interval as the
  scalar case. For the continuous construction the `gamma_k` are ordinary
  coefficients, the Jacobian column splitting into
  `-delta_k 1(x > psi_k) Z_j`, and a sub-term's penalty passes through:
  `seg(x, psi ~ random(~1 | id))` is the random-changepoint model of
  Muggeo, Atkins, Gallop and Dimidjian (2014). For `jump` the Fasola
  identity splits the `gW` column into `W Z_j`, whose coefficients are
  `c_k = -kappa_k gamma_k`, and the development is read off
  componentwise, `gamma_k = -c_k/kappa_k`, exactly as the scalar
  break-point is; a sub-term carrying a penalty is rejected there, the
  penalty acting on the development scaled by the jump size.

  `jseg()` REJECTS a development, and the reason is measured rather than
  presumed: its reading of the break-point is a quadratic in the
  increment that couples the slope change with the jump and does not
  split over the columns of a development, while the componentwise
  reading that remains diverges whenever the jump size passes near zero
  mid-iteration -- from the grid start on an ordinary sample the scalar
  quadratic settles at the truth and the linear reading runs to the
  confinement boundary. A development of `~1` reproduces the scalar
  construction: same block, same contribution, and the iteration walks
  the two to the same break-point (pinned to 1e-6).

* A subformula of `nl()` goes through `interpret_formula()`, so it takes
  any term of the package, and it is written as a two-sided formula in
  `...` whose left side names the parameter.

  `nl(~ a * exp(-r * x), a ~ ridge(~g))` is a population value (the
  intercept) plus shrunken departures; `a ~ s(z)` lets a parameter move
  smoothly with a covariate; `a ~ random(~1 | g)` is a random intercept on
  the parameter's unconstrained scale. The old route built the sub-design
  with `stats::model.matrix()` directly, so `a ~ g` worked and
  `a ~ ridge(~g)` could not exist. The penalties the sub-terms carry are
  reported through `term_penalties()` under the key `parameter::subterm`,
  with indices in the term's own numbering, so a fitting layer reaches
  their hyperparameters as it does any other term's. Because the
  development acts on the unconstrained scale, the parameter stays in its
  own set at every observation whatever the coefficients are.

  `subformulas = list(a = ~g)` remains as the programmatic spelling; a
  parameter may carry one subformula, whichever spelling supplies it. A
  structural term, and a term whose block moves with its own coefficients,
  are rejected: a parameter's submodel must be a fixed design. Prediction
  reapplies each sub-term's own blueprint (levels, knots, constants)
  instead of the hand-rolled terms/xlev record the old route kept, and the
  Jacobian is assembled block by block so a sparse sub-design -- a random
  intercept's indicators -- stays sparse through the block it becomes.

* The charts of a score-driven term's parameters are configurable, and the
  score loadings ride the log link by default.

  `gas(links = list(alpha1 = identity_link()))` overrides the default of any
  base parameter; a deviation cannot be named, being unconstrained by
  construction. The defaults are now: the level on the identity, every
  loading on the LOG link, the persistence on the rhobit chart of its
  partial autocorrelations. A positive loading responds in the direction of
  the score, which is the case the literature writes, and positivity is
  structural: a deviation (or, later, a submodel) moves the loading on the
  log scale, so no group can take a negative one -- which under the old
  identity default it could.

  The chart's own curvature now enters `term_curvature()`: the first and
  second derivatives of the level and the loadings in their coordinates are
  the links' rather than one and zero, verified against `numDeriv` at
  p and q up to 2 through the existing curvature tests, which read
  `term_links()` and exercised the new terms the moment the default moved.

* `term_start()` says where a term's own parameters start, on the
  unconstrained scale.

  The start belongs to the term because only the term knows what a
  coordinate of zero means on each of its charts. The base method returns
  zero everywhere, each link's natural point; `gas()` overrides it for the
  loadings, whose log chart puts a loading of ONE at zero -- a response
  strong enough to destabilize the recursion at ordinary curvatures -- and
  starts them at 0.1 on the parameter scale instead, whatever chart each
  one rides.

# modelterms7 0.26.0

* A sparse matrix handed to `ridge()`, `lasso()`, `enet()`, `scad()` or
  `mcp()` is no longer densified.

  `.penalized_spec()` called `as.matrix()` on any non-formula input, so a
  caller's `dgCMatrix` was densified at SPECIFICATION time, before anything
  had a chance to keep it. Measured on a 4000 x 60 indicator design at
  density 0.017: **0.050 MB became 1.920 MB**, a factor of 1/density, and
  the block stayed dense from there on. A penalized block is exactly where a
  sparse design turns up -- indicators over many levels are what a lasso is
  for -- so this was the one input the constructor had to preserve and the
  one it destroyed. `term_predict()` did the same at new data.

  Both keep the Matrix now, and a logical Matrix is carried to double rather
  than rejected, an indicator being the commonest sparse input. The
  infrastructure was already there and unused: `.is_block()` has accepted
  sparse since 0.24.0 and `check_term()` already reported it.

  End to end through `statmod()`, against the same design densified by hand,
  with the coefficients identical:

  | | sparse | dense | |
  |---|---|---|---|
  | ridge, n=20000 p=200 | 0.75 s | 4.73 s | 6.3x |
  | ridge, n=20000 p=1000 | 2.37 s | 90.86 s | **38.3x** |
  | lasso, n=20000 p=200 | 1.04 s | 3.39 s | 3.3x |
  | lasso, n=20000 p=1000 | 3.61 s | 59.22 s | **16.4x** |

  The coefficients agree to 2.4e-16 for the smooth branch and EXACTLY for
  the kinked one.

  ⚠️ The gap between the two branches is not noise and is worth reading:
  the kinked branch reaches a compiled coordinate descent that takes an
  `arma::mat`, so `statmodels7:::coord_fit()` still materializes the
  penalized block dense at that one boundary. That is the remaining
  densification in the chain, it is the whole of the difference between 38x
  and 16x, and closing it means a sparse path in the kernel -- which is also
  the natural algorithm, a coordinate update on a sparse column touching
  only its nonzeros.

# modelterms7 0.25.0

* `ridge()`, `lasso()`, `enet()`, `scad()` and `mcp()` take `standardize`.

  A hyperparameter is comparable across coordinates only where the
  coordinates share a scale, and without this a lasso penalizes a column
  measured in metres more than the same column measured in kilometres.
  `standardize = TRUE` divides each coefficient by the standard deviation of
  its own column through the penalty's DIAGONAL MAP, so the design is never
  rescaled: a sparse block stays sparse, `lambda` stays one number, and the
  coefficients come back on the scale the data arrived in with nothing to map
  back. Centring, which is what would destroy sparsity, is not needed, the fit
  being invariant to a translation of a penalized column where an intercept is
  free.

  The spread is computed from the BUILT block and frozen in the blueprint, so
  the same term standardizes identically in every equation of a distributional
  model and does not move with the working weights of a fit. A constant column
  takes `s_j = 1`. `print()` shows the values, a number that changes the
  meaning of `lambda` having to be legible.

  ⚠️ For SCAD and MCP this is not a rescaling of `lambda`, and the naive
  substitution is wrong by a wide margin. Measured against the published
  piecewise forms transcribed independently, at spreads from 0.5 to 3:
  `rho(s b)` differs from `rho(b; lambda*s, a)` by **11.2** for SCAD and from
  `rho(b; lambda*s, gamma/s)` by **4.39** for MCP. The exact relations are
  `s^2 rho(b; lambda/s, a)` and `s^2 rho(b; lambda/s, gamma)` -- an overall
  factor as well as both hyperparameters -- and SCAD is not a member of its
  own family at any parameters without that factor. The diagonal map expresses
  `rho(s b)` exactly (0), which is why no new arithmetic was needed; a test
  pins all of it, including the two wrong answers.

* `random()` does not standardize and takes no such argument, its columns
  being grouping indicators and its hyperparameter a variance component;
  `gas(deviations =)` likewise, a deviation being a parameter of the recursion
  rather than a coefficient on a column. Both reject the argument by name.

# modelterms7 0.24.0

* A grouping indicator is built SPARSE. `.random_block()` assembled a dense
  `n x (m*d)` block, and its intermediate `outer(g, levels(g), ==)` was a
  dense `n x m` besides. A row belongs to one group, so the density is `1/m`
  whatever the data. At `n = 20000` and `m = 1000`: **152.6 MB against 0.23
  MB**, built in 1.76 s against 0.0011 s, and the crossproduct every fitting
  iteration takes in **0.0006 s against 12.77 s**.

  `check_term()` asks `.is_block()` rather than `is.matrix()`, which is FALSE
  for every `Matrix` class and would have failed a term for being efficient.
  A term's `X` was already `class_any`, so the contract needed no change.

* `term_curvature()` accumulates PER GROUP and never forms the square over
  all the unknowns. A group's rows reach the coefficients, the population
  parameters and that group's own deviations, and nothing else, so the active
  set has the same size whether the panel has ten groups or a thousand.
  Measured against the full square, which it reproduces EXACTLY:

  | groups | unknowns | full square | per group | speedup |
  |---|---|---|---|---|
  | 25 | 79 | 0.31 s | 0.100 s | 3.1x |
  | 50 | 154 | 1.89 s | 0.060 s | 31.5x |
  | 100 | 304 | 12.48 s | 0.130 s | 96.0x |
  | 200 | 604 | 53.51 s | 0.390 s | 137.2x |

  `max|W_full - W_group|` is 0 at every size, not merely small: it is the
  same sum over the same terms, restricted where the rest is zero by
  construction. The cost per observation is now flat in the group count (37
  to 43 microseconds at 100, 250 and 500 groups), and what remains is the
  `n x m` jacobian, which is the return value and is O(n m) rather than the
  O(n m^2) the square cost.

  `blocks` is called with the row of the jacobian RESTRICTED to the active
  set and with that set, returning its pieces in the same coordinates. A
  three-argument callback of the earlier shape still works and is given the
  full row, paying the quadratic allocation the restriction avoids.

# modelterms7 0.23.0

* A score-driven term carries the names its literature uses. The score
  loadings are `alpha1`, `alpha2`, ... where they were `a1`, `a2`: they are
  the quantities themselves, each on the identity link, so the name can
  promise what it reports.

  The persistence is a different matter and keeps its own name. It rides a
  PARTIAL AUTOCORRELATION, the stationary region of an autoregression not
  being a box, so a free coordinate called `beta1` would promise the
  coefficient and report the chart -- and above `q = 1` the two are
  different numbers.

* `term_readable()` is the new generic that reports what a fitted term is
  about, with the Jacobian from the term's own parameters, in the shape
  `parameters7::param_readable()` already uses for a matrix parameter. The
  base method answers with the parameters themselves on the parameter
  scale, so every existing term is unchanged.

  A score-driven term answers with `omega`, the loadings, and the
  AUTOREGRESSIVE COEFFICIENTS `beta1`, `beta2`, ..., taken through the
  Levinson-Durbin recursion whose Jacobian the filter already computes and
  chained onto the rhobit link of each coordinate. Against `numDeriv` the
  Jacobian agrees to 1e-11 at every order tried; at `q = 1` the coefficient
  is exactly the link's inverse, which is the case where the two coincide.

# modelterms7 0.22.0

* A penalty is an object, not a string. `gas()`, `nl()`, `seg()`, `jump()`
  and `jseg()` took `"none"`, `"lasso"` or `"ridge"` through `match.arg()`,
  which put a term's reach at two of the penalties `penalties7` offers and
  made every other one need a name invented for it here. They now also take

  * a `penalties7` penalty, used as it stands, so an elastic net, a
    heavy-tailed prior or a structured precision reaches a term directly;
  * a function of the number of coefficients returning one, which is what a
    penalty whose WIDTH the data decide needs -- a panel's deviations exist
    only once the groups are counted, so a specification cannot name one.

  The two shorthands keep working and remain the defaults. A penalty given
  as an object is checked against the count where that count first exists,
  at `term_build()`, rather than being evaluated at a coefficient vector of
  another length and recycled in silence.

  `.penalty_factory()` is the one place that reads the argument and
  `.penalty_arg()` the one that validates it, so the three constructors
  cannot drift apart.

# modelterms7 0.21.0

* `term_curvature()` carries deviations, where it used to reject them.

  A group's parameters are the population values plus that group's
  deviations ON THE UNCONSTRAINED SCALE, which is the scale a deviation is
  defined on. That map is affine: its Jacobian is a matrix of ones and
  zeros and its second derivative is exactly zero. So the recursion is
  unchanged and only the lift widens, a base coordinate reaching two
  columns of the caller's unknowns instead of one, the population value
  and that group's own deviation. Nothing had to be derived.

  Against `numDeriv` on the filter itself, over deviations on every
  parameter, on the level alone and on a mixed pair, at `p` and `q` up to
  two: the Jacobian agrees to 3e-10 and the curvature to 5e-10 relative.
  Where the gap looks larger the reference is the weaker side -- at
  `q = 2` two Richardson settings disagree with each other by 4.4e-4 while
  the closed form sits 2.8e-7 from one of them, which is 2e-10 of the
  matrix's own size. Two structural claims are asserted as well: at a zero
  deviation the population block reproduces the shared-parameter term's
  curvature, and whatever the deviations are, one parameter's deviation
  columns sum to its population column, which is the affine lift.

  This is what a penalty over a panel's deviations needed: with it,
  `statmodels7` fits such a term, inverts its information and estimates
  its hyperparameters.

# modelterms7 0.20.0

* `term_hessian()` returns the exact Hessian of a likelihood mixed over
  latent states, in the whole of a caller's unknown vector: the
  coefficients of every equation together with the term's own parameters.

  What a caller could assemble from `term_posterior()` alone is the
  COMPLETE-DATA information, the ordinary one averaged over the smoothed
  states. That is the matrix an EM step inverts and it is not the observed
  information: by the missing-information principle the two differ by the
  conditional variance of the complete-data score, so the complete-data one
  is the larger and a standard error read off it is too small. Measured on
  a two-regime gaussian, the difference reaches 30 per cent of a standard
  error where the regimes overlap and vanishes as they separate and the
  states become known.

  Louis's identity is one route to it and is not the one taken. The scaled
  forward recursion computes the observed log-likelihood exactly, as a sum
  of the logarithms of its normalizing constants, so differentiating that
  arithmetic twice gives the observed Hessian with no identity, no pairwise
  smoothed probabilities and no second-moment recursion. The first and
  second derivatives of the filtered distribution are propagated beside it
  and renormalized by the quotient rule. Louis's identity becomes the
  check instead: the difference between the two matrices must be positive
  semidefinite, and it is measured to be, strictly, wherever the states
  carry any uncertainty.

  The cost is `O(n K^2 m^2)`, and the computation is meant to run once at a
  fitted point. `regime_stationary()` gained the second derivative of the
  stationary distribution, from the same linear system as its first.

# modelterms7 0.19.0

* `term_posterior()` returns the smoothed state probabilities of a latent
  Markov term, which is everything a model layer needs to differentiate a
  likelihood mixed over states. By Fisher's identity the derivative of that
  likelihood in ANY predictor the model carries is the posterior-weighted
  derivative of the ordinary one, so a caller differentiates its own
  log-density K times vectorized and weights, and needs no callback per
  observation. That is the property that made the forward pass compilable
  read once more: a regime shifts a predictor known before the recursion
  starts.

  The probabilities come from the forward pass this term already runs and a
  backward pass beside it, both normalized -- without which the quantities
  are products of t densities and reach zero in double precision within a
  few hundred observations. Validated against `numDeriv`: the rows sum to
  one to 1.1e-16 and Fisher's identity holds to 8.1e-09, 1.3e-07 and
  7.7e-10 over one series, three regimes and groups, where the score at the
  marginal mean is out by 1.4.

* `term_level_param()` says which of a term's parameters shifts its
  equation's predictor by a constant: `"omega"` for `gas()`, `"level1"` for
  `regime()`, and `character(0)` for everything else. It exists so that a
  fitting layer can resolve the confounding with an intercept rather than
  refuse the model. Which parameter is the level is the term's answer;
  which one is dropped is the layer's, since only the layer knows what else
  the equation carries.

# modelterms7 0.18.0

* `term_curvature()` is the second-order companion of `term_adjoint()`: the
  forward Jacobian of the predictor a structural term produces in a caller's
  unknowns, and the second derivative of that predictor contracted against
  the caller's weights. It is what an observed information needs and what
  the reverse recursion alone does not give.

  The contract keeps the split the adjoint already uses. `seed` is the
  derivative of the static predictor in the caller's unknowns, so the term
  learns nothing else about them; `blocks` is a callback returning, at an
  observation and the Jacobian the recursion has reached, the two model
  quantities that seed the first and second derivatives of the score --
  `sum_q l_pq C_q` and `sum_{r,r'} l_prr' V_r' V_r'`. A model of one
  equation supplies zero and `l_ppp D'D`.

  Against numDeriv at p, q in {1,2}^2: the Jacobian 1.05e-10 and the
  contracted second derivative 1.0e-09 on a matrix of scale 24.6.

  Deviations are refused rather than silently mishandled: the per-group
  chain adds a factor to every derivative and is not written.

* `gas_levinson2()` carries the second derivatives of the autoregressive
  coefficients in the partial autocorrelations, which the curvature needs
  because the persistence reaches the predictor through that map. The
  recursion is bilinear, so differentiating twice adds no new kind of term,
  only the two places the product rule puts the first derivative. Against
  numDeriv: 3.8e-12, 2.9e-12 and 9.4e-11 at q = 2, 3 and 4. The value and
  the jacobian are bit-identical to `gas_levinson()`, which is what says
  this is not a second route to them, and the last coefficient's second
  derivative is exactly zero at every order, it being the last partial
  autocorrelation.

# modelterms7 0.17.0

* `nl()`'s numerical route uses numericals7's stencil library. It wrote out a
  three-point central difference with its own step, and modelterms7 did not
  import numericals7 at all -- the nodes, the weights and the step are
  `fd_offsets()`, `fd_weights()` and `fd_step()` now, at accuracy four, which
  is the five-point rule. Measured against the exact symbolic Jacobian of the
  same function, on the route an opaque `f(x, theta)` takes: 3.7e-13 against
  the three-point rule's 2.7e-11. This Jacobian is the design block, so its
  accuracy is the accuracy of every step such a fit takes; the two extra
  evaluations are the trade distributions7 measured for the skew t.

* `term_adjoint()` is the reverse recursion of a structural term: the
  derivative of a caller's objective with respect to the static predictor the
  term was handed, and with respect to the sequence of scores it was given.

  `term_filter()` returns the derivative of the predictor in the term's OWN
  parameters, which is what estimating those needs, and it is not what
  estimating the coefficients of the same equation needs. A score-driven
  level at one time is driven by the scores at earlier ones, read at
  predictors those coefficients also enter, so the derivative of the
  predictor in a coefficient carries a term the block does not. Measured
  against `numDeriv` on the derivative in the static predictor: the reverse
  recursion agrees to 1e-8 and the direct score alone is wrong by 0.6 to 1.05
  in every configuration tried -- one series, p = 2 and q = 2, groups, and
  groups with deviations.

  Propagating that forward would cost one derivative array per coefficient;
  the reverse pass costs one whatever their number. Two derivatives are
  returned rather than one because the score depends on more than the
  predictor it is read at: multiplying `dscore` by the mixed second
  derivative of a log-density gives the derivative in ANOTHER equation's
  predictor, which is what a model layer with several distribution parameters
  needs.

* `term_value()` takes `newdata`, so a fitting layer can compute a term's
  contribution on other rows. Where the block is a Jacobian, `term_predict()`
  times the coefficients is the linearization and not the contribution: for
  `seg()` the two differ by a step at the break-point in a construction that
  is continuous. Rows are treated as `term_predict()` treats them, through
  the levels and constants the blueprint recorded.

* `term_converged()` asks whether a term's own iteration has settled, which a
  score cannot always answer. Where the block is the Jacobian of the
  contribution the gradient of the model's objective is the model's and its
  vanishing is the test; where the block is a working linearization with a
  frozen weight, as in `jump()` and `jseg()`, the profile objective is a step
  function in the break-point and has no gradient to vanish. The base method
  is `TRUE`, and the segmented method is `seg_converged()`.

# modelterms7 0.16.0

* `gas()` carries a population value and a deviation per group. With
  `by` and `deviations` each group of a panel is filtered with parameters
  of its own, written as

      psi[j, i] = g_j^-1( g_j(psi_j) + delta[j, i] ),

  the deviation acting on the unconstrained scale of the chart the
  parameter lives on, so a persistence stays inside (-1, 1) whatever the
  deviation is. `deviations` takes TRUE for every parameter or the names
  of the ones that carry one, and needs `by`. The deviations are
  parameters of the term, named after the parameter and the level
  (`omega.dev.a`), and carry the identity link, being unconstrained
  already.

  They are parameters and NOT a penalty on the per-group values through a
  difference matrix, which is the same model written the other way. The
  difference decides what can be fitted: a penalty over a general map is
  the generalized-lasso problem, whose proximal operator does not split by
  coordinate, while a deviation named as a coordinate is reached by a soft
  threshold and by a coordinate descent unchanged.

  The filter runs once per group and chains the columns of its jacobian
  onto the population values, exactly, `d psi[j,i] / d psi_j` being
  `g^-1'(g(psi_j) + delta) g'(psi_j)` and the derivative in the deviation
  the same without the second factor. At a zero deviation the two are
  reciprocal by the inverse function theorem, so the filter and every
  population column are then bit for bit the shared-parameter ones. The
  jacobian agrees with `numDeriv` to 1e-10 with deviations on one
  parameter and on all of them.

  The deviations are identified by their penalty and not otherwise. A
  parameter and its m deviations are m+1 numbers describing m group values,
  so a constant added to the population value on the unconstrained scale and
  subtracted from every deviation leaves the filter exactly unchanged, and
  the likelihood is flat along one direction per parameter carrying them.
  That is the parametrization of a random effect, identified there by a
  variance component and here by the penalty, which selects the deviations
  of smallest size among the descriptions of the same panel.

* `gas(penalty =)`, `nl(penalty =, penalize =)` and `seg(penalty =)`
  declare the parameters they penalize through `term_penalties()`, naming
  the coordinates each penalty covers instead of selecting them from the
  block with a map. The map was the defect: a separable penalty under a
  selection map is the generalized lasso, so `penalties7::has_prox()` was
  FALSE for a `seg(penalty = "lasso")` and neither a proximal step nor a
  coordinate descent could be taken on it. Named as coordinates the map is
  the identity and both are available.

  - `gas()` penalizes the deviations, one penalty per parameter carrying
    them, and rejects a penalty without them: the population parameters of
    a filter are not shrunk towards zero.
  - `nl()` penalizes one parameter at a time, the whole coefficient vector
    where the parameter carries a subformula, so a lasso there selects
    which covariates a parameter depends on. What is shrunk is the
    coefficient, so with a link the target is `g^-1(0)` and not zero.
  - `seg()` penalizes the changes as before, and `jseg()` now declares two
    penalties rather than one over their union: a slope change and a jump
    are not comparable quantities and cannot share a hyperparameter.

* `edf()` counts a term parameter by parameter rather than reading one
  penalty for the whole block. A parameter no penalty reaches counts one;
  a parameter under a kinked penalty counts one when it is away from zero;
  the rest are counted together by `tr[(H+S)^-1 H]` over the sub-block they
  occupy, with `S` carrying each smooth penalty's Hessian at the parameters
  it covers and zero elsewhere. Each rule reduces to what the term reported
  before when one penalty covers the whole block. `theta` is that penalty's
  hyperparameters for a term carrying one, and a list keyed by the penalty
  names for a term carrying several.

* `term_smooth()` reads `term_penalties()` too, so a term penalized over
  part of its parameters answers for the part: `seg(x, penalty = "lasso")`
  is not smooth although its linear effect and its break-points are
  unpenalized.

* `term_npar()` answers for a structural term as well, counting the
  entries of `term_params()`, which is the vector `term_penalties()`
  indexes into there.

# modelterms7 0.15.0

* `term_penalties()` is what a term declares it wants penalized: a list of
  entries, each naming a subset of the term's own parameters and the penalty
  over them. It generalizes `term_penalty()` in the two directions a model
  layer needs.

  A term may carry MORE THAN ONE penalty, over different parameters of its
  own -- a panel model with a population value free and a deviation per group
  shrunk is one penalty over part of the parameters and none over the rest.
  And the parameters need NOT be coefficients of a design block: the
  persistence of a score-driven term, the nonlinear parameters of `nl()`, the
  break-point of `seg()` are parameters of the term and nothing else, and all
  a penalty needs of them is a vector of numbers and their positions.

  The base method answers from `term_penalty()`, so every term shipped here
  needs no method of its own and behaves exactly as before, and a structural
  term answers with an empty list rather than raising, which is what lets a
  caller enumerate over every term without knowing which kind it has.

  The entry's name is unique WITHIN the term and is empty for a penalty over
  the whole of it. It is not the term's name: two `ridge()` terms in one
  formula are two terms with their own hyperparameters, and it is the caller
  that knows what it called each one.

# modelterms7 0.14.0

* te() centers its block: the tensor product of the marginal bases
  contains the constant, and the null space of the tensor penalty
  contains it too, so beside an intercept the design was rank deficient
  by exactly one with no penalty covering the deficiency. Measured on
  300 observations of te(a, b, k = 5), the design of y ~ te(a, b, k = 5)
  had 25 of 26 columns, a smallest singular value of 2.2e-15 and a
  condition number of 8.0e15, and the penalized information a smallest
  eigenvalue at the rounding floor -- which chol() accepts or rejects by
  the luck of rounding, so vcov(), confint() and the outer criterion
  returned numbers computed on a singular matrix.
  The block now carries the sum-to-zero constraint over the observed
  covariates, through basis7::constrain_basis() applied to the product
  basis, as mgcv does for a smooth. The same design is 25 of 25 columns
  at a condition number of 138.5, and the penalized information has a
  smallest eigenvalue of 0.93.
  A tensor term has one column fewer than the product of its marginal
  dimensions: te(x, z, k = 4) reports 15 parameters where it reported
  16. The penalty follows by congruence and its rank does not move, the
  direction removed having been one of its null directions (21 of 24
  where it was 21 of 25). The transform is stored in the blueprint and
  reapplied by term_predict(), as the Demmler-Reinsch transform of s()
  is. The level of the surface is the model's intercept, so a formula
  removing it fits a surface constrained to average zero.

# modelterms7 0.13.0

* enet(): the elastic-net term, beside ridge, lasso, scad and mcp,
  carrying penalties7::elasticnet_penalty(). Like the lasso it is not
  smooth, and its effective degrees of freedom are the nonzero count.

# modelterms7 0.12.0

* regime()'s forward recursion is compiled (src/regime_forward.cpp),
  and the density and score of every observation under every regime
  are computed once by k vectorized calls instead of 2nk scalar ones.
  Unlike the score-driven filter, nothing has to call back into R: a
  regime shifts the predictor by a level of its own, so none of those
  values depends on the filtered state. Measured against the R form,
  kept as the twin .regime_forward_r: 4.5x at k = 5, 13x at k = 3 and
  28x at k = 2, over T from 1e3 to 1e5, agreeing to 1.8e-15.
  End to end the term costs 2.04 microseconds per observation at
  k = 3, against 40 before.
* term_loglik()'s closures are called with the whole index vector.
  A closure returning one value where n were asked is rejected with a
  message saying so.

# modelterms7 0.11.1

* interpret_formula() rejects a call that evaluates to neither a model
  term nor a covariate, naming the call, its class and -- when the
  function that was called is not the one modelterms7 exports under
  that name -- the package that masked it. mgcv exports s() and te()
  and segmented exports seg(), so a user with either attached wrote
  our formula and got theirs; the value used to travel to
  model.matrix and fail there, naming neither the call nor the mask.

# modelterms7 0.11.0

* The discontinuous terms follow Fasola, Muggeo and Kuchenhoff's
  algorithm as published. The weight is no longer capped: the
  covariate is rescaled away from the break-point by a factor `c0`,
  which leaves a gap around it, and `term_refresh()` halves the
  factor whenever the break-point reverses direction. `band` is
  gone and no damping is needed, the factor governing the step.
  Measured over eight samples, `jump()` now recovers the
  break-point from every starting position tried, where the capped
  form recovered it only from within a narrow basin.
* seg_step() and seg_converged() report the progress of the
  iteration and Fasola et al.'s stopping rule.
* seg_start() chooses the starting positions by scoring an equally
  spaced grid on the least-squares profile, which is the
  initialization Fasola et al. recommend and the piece that was
  missing. It is what settles the joint term: measured over eight
  samples of a jump and a change of slope at the same point, a
  single conventional start recovers the break-point in none to
  half of them depending on where it is placed, and the grid in all
  of them. Bootstrap restarting was measured beside it and does far
  less (0.12 to 0.75).
* Measured against the segmented package, on one covariate with one
  and two break-points and n from 200 to 20000: the continuous case
  agrees to four significant figures on the residual sum of squares
  and runs 2.1 to 5.5 times faster in 4 to 7 iterations; the
  discontinuous case agrees where both start inside the basin and
  runs 1.1 to 1.8 times faster.

# modelterms7 0.10.0

* seg(), jump() and jseg(): break-points estimated with everything
  else. The continuous case is the Jacobian of Muggeo (2003); the
  discontinuous one uses the identity of Fasola, Muggeo and
  Kuchenhoff (2018), which is linear in the break-point once the
  weight is frozen, so the break-point is read off two coefficients
  rather than searched for. Both run under the existing
  term_refresh() contract, with by and an optional penalty on the
  changes. The working block is compiled (1.2x to 3.2x the R form
  over n from 1e3 to 1e6, agreeing to a rounding); the linear fit
  around it is BLAS in either language.

# modelterms7 0.9.0

* regime(): a latent Markov chain of regimes, each shifting the
  predictor by a level of its own, with the likelihood evaluated by
  the normalized forward recursion and its derivative propagated
  beside the state. Built on parameters7::transition_matrix(), which
  had no consumer until now.
* term_loglik(): the second shape of the structural branch, for a
  term whose contribution is a likelihood rather than a predictor.

# modelterms7 0.8.0

* nl(): the nonlinear parametric term. The design block is the
  Jacobian of the contribution in its parameters, refreshed by
  term_refresh() as they move, and term_value() reports the
  contribution a Gauss-Newton step needs beside it. The function may
  be a formula, read symbolically where deriv() manages it and
  differenced where it does not, or an opaque function, always
  differenced. Links per parameter; covariate submodels on the
  formula route, which is the only one that says where a parameter
  enters.

# modelterms7 0.7.0

* The score-driven recursion is compiled. The two callbacks into R
  remain, the score and the curvature belonging to the model's
  distribution, but the arithmetic around them was 73 to 83 per cent
  of the loop's time: measured 2.3x to 3.2x faster, and the R loop
  stays as the twin the kernel is compared against.

# modelterms7 0.6.0

* The structural branch is real: term_params(), term_links() and
  term_filter() define what a term that rewrites the likelihood must
  provide, and gas() implements score-driven dynamics over groups and
  time, with the persistence carried on a partial-autocorrelation
  chart and the derivative of the filter propagated alongside its
  state.
* s() and te(): penalized smooths of one and several covariates, the
  first under the Demmler-Reinsch reparametrization that separates the
  linear effect from the nonlinear deviation, both accepting `by`.

# modelterms7 0.5.0

* random() accepts slopes (~ x | g) beside intercepts, with the
  within-group Gaussian unstructured or diagonal by `correlated`, a
  per-group precision structure replicated across groups through
  parameters7::kron_identity(), or a distribution applied
  coordinatewise.

# modelterms7 0.4.0

* random(~ 1 | g): grouped random intercepts with the effect
  distribution as the penalty -- independent Gaussian by default, a
  parameters7 precision structure, or a distributions7 object applied
  coordinatewise. Random slopes are rejected pending the
  block-diagonal composition in parameters7.

# modelterms7 0.3.0

* edf() with the counting rule per penalty (exact count, trace of the
  penalized smoother block, nonzero count), the penalty shown by
  print() on a built penalized term, and plot() at fitted
  coefficients.

# modelterms7 0.2.0

* The penalized terms ridge(), lasso(), scad() and mcp(), over formula
  or matrix input, with the penalty attached at build time and the
  smoothness flag read from its kink set.

# modelterms7 0.1.0

* First release: the term classes and generics, the unpenalized
  parametric term, the formula interpreter with recognition by
  evaluation, the censored-response constructor, and the term validator.
