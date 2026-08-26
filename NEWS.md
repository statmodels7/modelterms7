# modelterms7 0.60.0

* `term_simulate()` says how a response is DRAWN from a term that carries
  state, which is a different operation from fitting one and differs in
  which direction the response moves. A term whose contribution does not
  read the response reports it and leaves the drawing to the caller; a
  score-driven term cannot, its level being driven by the score of the
  response at the time before, so it draws as the recursion runs. One
  contract covers both: the caller passes a function that draws at a
  predictor, and a method that drew returns the responses while one that did
  not returns `NULL`. The score-driven method writes NO new recursion --
  `term_filter()` calls its score callback exactly once per observation, in
  time order, at the predictor just produced, so a callback that draws there
  turns the filter into a generator. A latent chain draws its path from the
  STATIONARY law the likelihood is written with, any other start being a
  different model from the one a fit reads back; a marginal break-point term
  draws each group's positions from their prior.

* `term_continue()` says what a structural term's contribution does at rows
  that come after the ones it was built on. A term with state cannot be
  reapplied -- what it reports at one row is where a recursion has got to --
  so a prediction past the series carries the state forward instead. What
  makes that possible without simulation is the model's own defining
  property: the quantity driving the recursion has zero conditional mean, so
  beyond the data the recursion is deterministic, `f` decaying towards
  `omega / (1 - beta)`. A row is placed by its own time within its own
  group; one falling inside the observed series is rejected with the rows
  named, that being a re-reading rather than a continuation, and so is a
  group the fit never saw. Validated against the filter re-run over the
  extended series, which shares the recursion and not the continuation:
  identical to the bit, alone, beside a covariate, over a panel and with the
  level developed over covariates. The base method signals an error rather
  than returning zero, which would read as a term with no effect.

* `term_static_deriv()` says how the predictor a structural term produces
  moves when the static part of the predictor moves. A score-driven level is
  driven by scores read AT the predictor the recursion is producing, so a
  coefficient in the same equation reaches the level as well as the design
  row; the derivative obeys the recursion the filter already runs, seeded at
  zero because the starting level depends on the term's parameters alone,
  and it needs no callback, the curvature it multiplies being the one
  `term_filter()` returns. Measured on a score-driven mean with one
  covariate beside it, a standard error that counts the static row alone is
  about a quarter too small; with the propagation it agrees with a numerical
  derivative of the predictor to 1.8e-09. The base method returns `NULL`, a
  term without state carrying nothing to propagate.

* `term_components()` says how a term's columns divide among the parameters
  it is written in: one entry per parameter, with the columns that belong to
  it and the sub-terms developing it. What divides them is the term's answer
  and cannot be recovered from a coefficient name, which is built for a
  reader; the base method returns an empty list, the honest answer for a term
  whose columns are one block with one meaning.

* Each entry also carries `sub_index`, the columns belonging to each of its
  sub-terms. A parameter may be developed by several at once and they need
  not be of one kind -- `seg(x, psi ~ random(~1 | id))` develops the
  break-point with an unpenalized intercept AND a random block -- so a
  consumer that reports a component reports a sequence. The division rests on
  how the block is assembled, which is the term's business, and a consumer
  that computed it from coefficient counts would be assuming it.

* `term_readable()` refuses quantity by quantity rather than all at once. A
  developed parameter is a vector of coefficients over covariates and has no
  single value to report, so it is skipped; the parameters beside it are
  unaffected. Refusing them along with it left a summary printing the
  working coefficients a discontinuous construction is fitted through, which
  are no part of the model. Where every parameter is developed nothing is
  left and the answer is `NULL`, as before.

* `term_components()` answers for a score-driven term. A structural term
  contributes no design columns, so what divides there is the PARAMETER
  vector: `index` gives positions in `term_params()`, which is the vector its
  state, its readable quantities and its variance matrix are all indexed by.
  In both cases the field names the term's own coefficients.

# modelterms7 0.59.1

* The zeros the sharp break-point terms answer are pinned by tests rather
  than only argued for. For `seg()` the second derivative is exactly zero --
  not small, zero -- on every observation a perturbation does not carry
  across the break-point, at each of three steps. The one exception is a
  point mass and is tested as one: with an observation sitting exactly on
  the break-point the difference on that row grows fourfold as the step is
  quartered, which is what "almost everywhere" sets aside. Without that
  second half the first is satisfied by a sample with no observation near
  the break-point, which asks nothing.

* And the two kinds of zero are told apart. For `seg()` the block's own
  difference converges onto `term_block_deriv()` to 1e-10, so its first
  derivative exists and the second is what vanishes; for `jump()` and
  `jseg()` `term_block_deriv()` answers zeros while the block moves by five
  orders of magnitude at the scaling factor the annealing descends to,
  because there the block is a working linearization with a frozen weight.
  The zeros are a refusal in one case and a description in the other.

# modelterms7 0.59.0

* The break-point terms answer `term_block_deriv2()`. With `smoothed` an
  `abs_smoother` the block is the true Jacobian, so the closed forms are the
  smoother's own one order further up than `term_block_deriv()` reads them,
  and the only new quantity is the smoother's fourth derivative: the block
  reads it to order two, the first derivative to order three, this to order
  four. `abs_smoother` carries those orders as functions, so nothing in
  `penalties7` had to change.

* The sharp constructions answer zeros, and the two reasons are different.
  For `seg()` the second derivative really is zero away from the
  break-points, the truncated line's derivative in the position being an
  indicator whose own derivative is a point mass, and the position column
  being linear in the change. For `jump()` and `jseg()` the block is a
  working linearization with a frozen weight rather than a Jacobian, which
  is why the first-order generics already answer zeros there.

* Two exact properties of the smoothed branch are pinned by tests. Where a
  break-point sits against its confinement limit the whole contribution is
  zero, every addend carrying a direction in the break-point -- which the
  first derivative does not, its position column moving with the change
  whatever the position does. And under `smooth_quintic()`, exact outside
  the width, the answer is zero on every observation further than the width
  from a break-point; that smoother is C^3, so its fourth derivative jumps
  at the width and the answer is exact away from those two points rather
  than everywhere.

* `.seg_smooth_parts()` takes an `order`, so the smoother's fourth
  derivative is computed only where this method asks for it and not at every
  step a fit takes.

* Still nothing consumes the generic. Measured against the previous release
  on three shapes -- a moving block, a fixed design, and a smoothed jseg
  with a penalized subformula -- the coefficients, log-likelihood, effective
  degrees of freedom, hyperparameters, criterion, outer gradient, outer
  Hessian and variance matrix are identical to the bit.

# modelterms7 0.58.0

* `term_block_deriv2()`, the second derivative of a design block in its own
  coefficients, contracted in two directions the caller supplies rather
  than returned as an array. It stands one order above
  `term_block_deriv()`, and the base method answers zeros, which is exact
  for a design that does not move with its coefficients and covers
  `linpar()`, `s()`, `te()`, `random()` and the five penalized constructors
  without a method of their own.

* `nl()` implements it. Nothing new is derived: the third derivative of `f`
  in its parameters is the third order of `nl_fderiv()` and the third
  derivative of the inverse link is `linkfunctions7::d3linkinv()`, so the
  method is an assembly of five addends, one for each way a factor of
  `term_block_deriv()`'s own expression can be differentiated again. A
  subformula is carried by the same expression, its design entering every
  addend, and a sparse development is scaled in its own storage rather than
  densified first.

* Nothing consumes it yet. The quantity enters the HESSIAN of a marginal
  criterion and nothing else, so no fit, no criterion value and no
  criterion gradient changes; what will change when `statmodels7` reads it
  is the standard error of a hyperparameter on a model whose block is
  curved, and the Newton direction of an outer search.

* `.nl_theta()` takes an `order`, so the third derivative of the link is
  computed only where this method asks for it and not on every Jacobian
  evaluation of every fit.

# modelterms7 0.57.0

* The two filter kernels and the curvature kernel take the shape
  `distributions7`'s `d7::par_for()` settled on. The worker's loop is
  noinline and the sequential branch runs through the worker rather than
  through a loop of its own, so both branches execute one compiled copy
  instead of two the compiler may optimize apart; and the worker installs
  the calling thread's floating-point environment before its chunk, which
  is what keeps R's `psigamma`, `bessel_k`, `pgamma` and `pbeta` from
  returning per-thread last bits. A filter accumulates over time, so a
  difference of that kind would be carried forward by the recursion rather
  than confined to the element that produced it.

* `threads` is passed to `parallelFor()` instead of being left to
  `RCPP_PARALLEL_NUM_THREADS`. `term_filter(threads = 2)` called outside a
  fit was running on every core of the machine; a fit, which sizes the pool
  through `numericals7::local_threads()`, is unaffected.

# modelterms7 0.56.1

* The gas fast-route twins compare the C route against the R callbacks at a
  tolerance of 1e-13 instead of `identical()`: clang on the arm64 macOS
  runner contracts the scalar composition's multiply-adds into FMAs, which
  R's interpreter never does, and the last bit moved. The tolerance still
  fails a wrong composition by many orders; everything computed by one
  route both times -- the inert-context fallbacks, the threads over
  groups -- stays `identical()`.

# modelterms7 0.56.0

* `seg()`, `jump()` and `jseg()` take `marginal = FALSE` (the default: the
  construction exactly as before, to the bit) or `marginal = TRUE`, the
  marginal break-point terms of `piano_marginal.txt`:
  `kind(x, psi ~ random(~1 | g), marginal = TRUE)` treats each break-point
  as a latent variable per group and integrates it out of the likelihood.
  The prior is part of the likelihood -- its parameters are ordinary ones
  estimated by maximum likelihood, `term_penalties()` declares nothing --
  and the term is structural of the likelihood shape (`MarginalBreakTerm`,
  the contract of `regime()`).

* For `jump()` the integral is EXACT, and it is not taken over the cells:
  the side process S_t = {k : psi_k <= x_(t)} over a group's sorted
  observations is monotone on the subset lattice with independent
  coordinates, so it is a HIDDEN MARKOV CHAIN on the 2^K side patterns
  whose transition factors over the coordinates -- a flip is weighted by
  its interval's prior mass, a coordinate that never flips contributes its
  tail mass through the final state's survival product -- and the forward
  recursion costs n K 2^K where the cell sum costs (n+1)^K. Measured end
  to end at the same answers to the printed digit: K = 2 goes 13.2 s to
  2.5 s, K = 3 goes 721 s to 8.8 s, and K = 5, unreachable before, fits in
  94 s with every position recovered within 0.05 of the truth. The cap is
  eight break-points, priced by the 2^K components a fitting layer
  evaluates the family at, not by the recursion. With one break-point the
  prior may be any continuous distributions7 family with its location
  fixed at zero (`random(distrib = fixed(student_t1_distrib(), mu = 0))`),
  the masses and their derivatives riding the cdf surface built for the
  censored likelihoods. The jacobian rides the same recursion; the
  posterior side patterns and the flip-interval posteriors of
  `term_latent()` come from the matching forward-backward pair; and
  `term_hessian()` propagates first and second derivatives through the
  recursion itself, so the observed Hessian is analytic for every K and
  every prior -- the toolkit's own finite differences left the marginal
  terms entirely, the continuous kinds' prior rows closing in the affine
  node motion and the step kind's in the propagated chain. The
  one-break-point gaussian interval-sum route survives as the twin the
  tests hold the propagation to.

* For `seg()` and `jseg()` (one break-point, gaussian prior) the
  conditional is smooth within an interval and the integral runs on a
  fixed Gauss-Kronrod panel per interval
  (`numericals7::gauss_kronrod15()`), interior nodes fixed points of the
  data so the prior-parameter derivatives read the prior alone; the region
  below the data, where the hinge keeps moving, is covered by panels that
  follow the prior's bulk, whose node motion the jacobian carries, and the
  region above it contributes its closed tail mass. Measured on the real
  integrand, GK15 per interval reaches 3e-9 where a 4-node Gauss-Legendre
  rule is off by 0.2 to 0.95 on the log-likelihood.

* The contract: `term_loglik()` returns one-step predictive contributions
  with the exact jacobian (numDeriv 1e-10 on the step kind, 2e-5 on the
  quadrature kinds' scale), `term_posterior()` the component weights
  Fisher's identity takes (side patterns for the step kind, node weights
  for the continuous ones), `term_hessian()` the exact observed Hessian --
  the one-break-point gaussian step differentiates the interval sum twice
  and doubles as the control; every other configuration assembles the
  component blocks exactly and takes the prior's rows from one central
  stencil on the analytic full gradient, the license the toolkit's
  non-closed derivatives run on -- and `term_start()` reads a two-stage
  exact profile off the target: pooled positions and coefficients with
  per-group intercepts, greedy over a quantile grid with every local
  minimum carried through a per-group one-position refinement and the
  winner chosen on the final per-group-position fit. A per-group profile
  with the full design was measured overfitting its own noise (a Poisson
  panel's linear effect at -4.7 against a truth of 0.15), and the pooled
  single-position profile alone picks the wrong basin (rss 103.6 at
  c = 8.3 against 104.8 at the truth's 4.4).

* New generics `term_levels()` (the shifts of a likelihood-shaped term's
  mixture components -- a vector for `regime()` and the step kind, a
  per-observation matrix for the quadrature kinds) and `term_latent()`
  (posterior means and standard deviations of the latent break-points,
  the truncated-prior moments through `numericals7::mills_ratio()` for
  the gaussian and `truncated()` + `expectation()` for an explicit prior,
  a moment the engine cannot deliver reported NA -- a heavy-tailed prior's
  edge intervals keep no mean below one degree of freedom and no variance
  below two).

* `smoothed`, `c0` and `n_boot` are ignored with a message under the
  marginal; a mixed fixed/random set of break-points, `by`, developments
  of other coefficients, several break-points for the continuous kinds
  and a non-gaussian prior beyond the single-break-point jump are
  rejected with the reason.

# modelterms7 0.55.0

* `seg()`, `jump()` and `jseg()` take `smoothed = NULL` (the default: the
  construction exactly as before, to the bit) or a `penalties7`
  `abs_smoother`. With a smoother the step and the hinge are replaced by
  their smooth versions, every break-point becomes an ordinary parameter of
  a `C^infinity` model, and the block is the TRUE JACOBIAN:
  `term_jacobian_block()` answers `TRUE`, the term is fitted by
  Gauss-Newton like `nl()`, and there is no working parametrization, no
  auxiliary `g` coefficient and no scaling schedule (`c0` is ignored, with
  a message). A development of a break-point -- `psi ~ random(~1 | id)`, a
  penalized one included -- is then legal for every kind, the read-off that
  constrained the discontinuous constructions having gone.

* The smoother's width is resolved at build from the covariate's spacing
  (the median gap between distinct values, within groups where a
  break-point development supplies a partition; `per_group = TRUE` keeps
  one width per group) unless the object carries one, is checked against
  the derived floor, and is reported by `print()`: it is the width of the
  transition, the bent-cable reading.

* `term_block_contract()` and `term_block_deriv()` have closed forms for a
  smoothed term of any kind, one order up in the smoother's own
  derivatives; `seg_psi()`, `seg_relocate()`, `seg_polish()` (whose exact
  profile does not depend on the mollifier), the lineage relabeling and
  `n_boot` all carry over. The restarts remain necessary: smoothing rounds
  the local optima, it does not remove them, and a smoothed fit from a bad
  start has been measured converging to an absurd local optimum while
  reporting success.

# modelterms7 0.54.0

* New generic `term_jacobian_block()`: whether a term's design block is the
  exact derivative of its contribution in its own coefficients (`TRUE`, the
  base method and `seg()`'s case) or a working linearization with a frozen
  weight (`jump()` and `jseg()`). It is what a fitting layer routes on: a
  Jacobian block licenses a Gauss-Newton step inside the model's own
  objective, a frozen one belongs to the fixed-point iteration of Fasola,
  Muggeo and Kuchenhoff (2018), which is not a descent method on that
  objective and stalls under a sufficient-decrease line search.

* `term_refresh()` on a break-point term relabels crossed lineages: when
  the implied positions come out of order, the (change, level, position)
  triples are permuted so the positions are ascending, each carrying its
  own scaling factor and direction -- what `segmented` does at every
  iteration. Two break-points that cross have exchanged roles, and left
  alone they chase the same feature until the block collapses onto a
  collision (measured on three break-points: psi = (-0.67, 0.50, -0.48)
  from a poor start). The contribution is a sum over the break-points, so
  relabeling moves no value; the coefficients the term stores are the
  relabeled ones, and a caller continues from those. A term whose
  per-break-point coefficients carry a development is left alone.

* `seg()`, `jump()` and `jseg()` take `n_boot`, how many bootstrap
  restarts (Wood 2001) the fitting layer runs after the iteration first
  converges -- `segmented`'s own default device, and the default here is
  its 10. The term declares the value; running the restarts belongs to the
  fitting layer, as with a penalty's hyperparameters. 0 disables.

* The operations a restarting loop runs on: `seg_reheat()` puts the
  scaling schedule back at `c0` with the directions and step record
  cleared -- the schedule only tightens, so an iteration resumed from a
  converged fit inherits factors at their floor and cannot travel
  (measured: restarts without the reset returned the incumbent unchanged
  ten times out of ten); `seg_relocate()` places the break-points at
  given positions with the changes kept and the block rebuilt, the
  read-off returning exactly the positions given; `seg_polish()` is
  coordinate descent over the positions on the exact profile, least
  squares at fixed positions being an ordinary linear model, optionally
  weighted so a bootstrap resample's profile is swept the same way; and
  `seg_profile_rss()` reads that profile at the current positions, which
  is what lets two configurations be compared at the cost of two linear
  fits rather than two model fits.

# modelterms7 0.53.0

* `term_adjoint()` on a gas term takes `fast` and `threads` and passes
  them to the forward pass it re-runs; and both filter kernels now return
  `curv`, the curvature read at each predictor -- computed anyway for the
  jacobian -- so the reverse pass looks the sequence up instead of
  evaluating the callback a second time at the same points. With a covered
  context the adjoint evaluates no R callback at all; without one the
  callbacks run once per observation instead of twice. Results are
  identical to the bit either way, and the twin test proves the fast route
  with a callback counter.

# modelterms7 0.52.0

* The second-order recursion of the SUBMODEL route runs compiled
  (`src/gas_curvature.cpp`). What made it portable is the same asymmetry
  as the filter's: a curvature runs at a point the caller has already
  fitted, so the score and the curvature of the density arrive as LOOKUPS
  (`score_values`, `curvature_values`) and the layer's blocks callback as
  DATA (`blocks_data`: the mixed second derivatives, the third
  derivatives one column per parameter pair, the static jacobian rows and
  the filter's parameter index), and the loop touches no R API. Coverage
  is stated rather than blurred: second order only and constant
  autoregressive charts -- a direction (the third order) or a developed
  partial autocorrelation keeps the R recursion, which stays the
  reference and is compared as a twin at a tolerance (the seg_block
  rule). The GROUPS run over threads (`threads`), each group's
  contribution merged on the main thread in group order, so the result is
  identical across thread counts to the bit; a test asserts it, and a
  callback counter proves the compiled route is the one taken.
  `.gas_curv_prep()` scatters every per-observation quantity once per
  group with the R route's own expressions, so the prep is shared
  arithmetic and not a second derivation.

# modelterms7 0.51.0

* The score-driven filter's kernels take a FAST context and a thread
  count. With `fast`, both recursions (the scalar route and the general
  submodel route) resolve the scalar C entry points of
  \pkg{distributions7} and \pkg{linkfunctions7} once per call with
  `R_GetCCallable` and read the score and the curvature without touching
  the R API -- mirroring `distrib_kernel()`'s composition expression by
  expression, held bit-identical to the callback route by a twin test; a
  family or link the registries do not cover leaves the context inert and
  the callbacks run as before. With no R in the loop the groups run over
  THREADS (`threads`, from the caller's `n_threads()`): a group's filter
  is independent and writes its own rows, so no reduction is split and
  the result does not depend on the count, bit for bit, which the tests
  assert with `identical()` at 1 and 2 threads. `term_filter()` passes
  both through its dots, so the contract is unchanged.

# modelterms7 0.50.0

* **`term_coef_start()` takes a `target`**, the response on the scale of the
  predictor the term contributes to. It is the one thing a term cannot work
  out for itself: the term knows its formula and its charts, the fitting layer
  knows the distribution, the link and the equation. Every method already had
  `...`, so the argument is compatible, and a term with no use for it ignores
  it.

* **`nl()` estimates its own parameters from the data.** Where a `target` is
  given, the parameters the caller did not pin with `start` come from a least
  squares fit rather than from zero, which for a nonlinear term is a
  degenerate point and not a neutral one: it linearizes where the function was
  never meant to be evaluated. On a logistic growth curve the start goes from
  `phi = 1, theta = 0, sigma = 1` to `50.71, 10.08, 2.071` against a truth of
  50, 10, 2, in 0.7 s.

  Two things make it work. The grid is DETERMINISTIC, so a start does not
  depend on the caller's random seed; a Latin hypercube of the same size, tried
  first, put the scale anywhere between 4.9e-06 and 2.07 over twenty seeds. And
  the parameters the function is jointly AFFINE in are separated out and solved
  by least squares at each point rather than searched over, read off
  `stats::D` as the fixed point of "the derivative in the parameter names no
  member of the set". Over five shapes and fifteen samples each, that is 15/15
  everywhere; on a four-parameter logistic it is the difference between 13/15
  and 15/15, with the worst relative error falling from 6.36 to 0.076. An
  opaque `fn` gets no separation, since a function does not say where its
  parameters enter, and still recovers the truth.

* A `start` reaches a parameter's submodel through a least-squares solve
  rather than being written at its first column, which was right only where
  the development carries an intercept. Where it is coded full rank -- which is
  what `phi ~ 0 + lasso(~g)` gives -- the old route gave the starting value to
  one level and zero to the others.

# modelterms7 0.49.0

* `term_block_deriv(term, coef, v)` is the ADJOINT of
  `term_block_contract()`: the block's own derivative taken along a direction,
  one entry per observation and column, where the other contracts over the
  observations and answers per coefficient. Neither computes the other. The
  GRADIENT of a marginal criterion needs the contraction; its HESSIAN needs the
  direction, `dK/dbeta` being wanted there in the direction the mode moves
  rather than traced.

* It reads the same closed form and needs no derivative the other did not: for
  `nl` it is the SECOND derivative of `f`, not the third. `seg` has it on the
  same two pieces its contraction carries; `jump`, `jseg` and every fixed block
  answer zeros.

* Verified against the block differenced along `v` -- 3e-10 to 2e-11 against
  scales of 1.2 to 6.8 -- and against `term_block_contract()` through the
  adjoint identity `sum_ij A_ij (dX.v)_ij == v' contract(A)`, which shares no
  code with either and comes back EXACTLY zero.

* Nothing in \pkg{statmodels7} consumes it yet: correcting the outer Hessian
  needs the same piece in three places at once, and doing one of them made a
  weakly identified fit worse rather than better.

# modelterms7 0.48.0

* `seg()` implements `term_block_contract()`, where it inherited the base
  method's zeros. Those are right for a fixed design and wrong for a block
  that moves, and a consumer differentiating anything built from `X` needs
  them: measured through \pkg{statmodels7}, the outer gradient of a penalized
  `seg` was `6.6e-04` relative against a central difference of the criterion
  where an `nl` sits at `4.3e-10`, and is `8.2e-10` now.

* What is written out is what is BOUNDED. The truncated line `(x - psi)_+` has
  derivative `-1(x > psi)` in the break-point, and the break-point column
  `-gamma(x) 1(x > psi)` has that same indicator as its derivative in the
  CHANGE and zero almost everywhere in the break-point, the indicator being a
  step. A caller could not take the difference instead: measured at h, h/4 and
  h/16 the quotient reads 3.6e4, 1.4e5 and 5.8e5.

* ⚠️ `jump()` and `jseg()` keep the zeros, and the reason is their
  construction rather than the work being unfinished. Their position is READ
  OFF a product of the unknowns, `psi = -g/delta`, so a column's derivative
  runs through that read-off rather than through a development's design, and
  the weight `W = 1/(2|x~ - psi|)` they carry has an unbounded derivative in
  the break-point. Their block is a working linearization with a frozen weight
  rather than a Jacobian. A test pins the zeros so a partial implementation
  cannot arrive unnoticed.

* ⚠️ The reference is a brute-force `dX/dbeta` and NOT the criterion. At two
  break-points the criterion's own central difference disagrees with ITSELF by
  7.37 relative, so it reported the correct contraction as wrong by 0.92; the
  brute force puts all four spellings -- one break-point, two, a developed
  change and a developed break-point -- at `1e-9` to `1e-10` against scales of
  4 to 25.

# modelterms7 0.47.0

* `sparse` defaults to `NULL` in `linpar()`, `s()`, `te()` and the five
  penalized constructors, and the storage is then SETTLED AT BUILD from the
  size of the design rather than asked for. `TRUE` and `FALSE` override it,
  and an explicit `TRUE` is still refused where there is nothing to build on.

* **The threshold is the cells of the dense indicator part, not a count of
  levels.** The dense form holds `n` times its column count against one
  non-zero per row, so that product is what separates the two routes.
  Measured end to end on `y ~ 0 + g + s(x)` over eighteen combinations of
  sample size and level count, they cross at about `1e5` cells, and the rule
  predicts every one of the eighteen: at `n = 1000` the sparse route loses at
  every level count up to sixty (0.93x at `6e4` cells), at `n = 5000` it
  crosses between fifteen and twenty-five levels, and at `n = 20000` between
  six and ten (1.00x at `1.2e5`). The same threshold accounts for the large
  cases, four hundred levels at `n = 20000` being `8e6` cells and 43.75x, and
  for the negative one: a design carrying no factor has no indicator part, and
  forcing the storage there measures 0.66x to 0.90x.

* The SETTLED storage is what the blueprint carries, so `term_predict()`
  reproduces the build's kind rather than deciding again on however many rows
  the new data has.

* `random()` is untouched, its block being sparse by construction.

# modelterms7 0.46.0

* `nl()` takes `gradient`, `hessian`, `deriv3` and `deriv4`: functions
  `function(theta, data)` returning the derivatives of the nonlinear function
  in its own parameters, each a NAMED LIST keyed as \pkg{distributions7} keys
  its own derivative surfaces. They are independent, so the orders worth
  writing out by hand can be written and the rest left alone, and the route is
  chosen ONE ORDER AT A TIME: a supplied function, then symbolic
  differentiation of the expression, then one stencil applied to the highest
  order that IS analytic.

* **Writing the Hessian pays twice.** The third and fourth orders are then one
  difference away from an exact second rather than from the function. Measured
  on `a * exp(-r * x)` given as an OPAQUE function, so that nothing is
  symbolic, against the closed forms: with nothing supplied the orders come
  back at 8.4e-03, 4.62 and 1.96e+03; with the gradient and the Hessian
  written out, 0, 2.2e-12 and 8.9e-11.

* The component names are normalized here, so `r_a` and `a_r` are one
  component and the order they are returned in does not matter. ⚠️ The
  accepted spellings are BUILT by permuting the parameter order and never
  obtained by splitting the user's string: a parameter whose own name contains
  an underscore makes a name ambiguous to read back, which is the trap this
  toolkit records for Hessian component names. A name that is not a component,
  a missing one or a repeated one is an ERROR at `term_build()`, an exact
  derivative that is quietly not used being worse than none.

* `nl_fderiv(term, coef, order)` reads any order back, and
  `term_block_contract(term, coef, A)` is the new generic a fitting layer
  needs: the contraction of the block's own derivative,
  `sum_ij A_ij dX_ij/dbeta_c`, in closed form at O(nm) for `nl()` and zeros on
  the base class, which is exactly right for a fixed design. The chain rule
  onto the coefficients -- each parameter's link and its subformula's design --
  stays in the term, the only thing that knows them. Verified against a
  brute-force numerical `dX/dbeta` at 1e-11 across links on one or both
  parameters, a developed parameter, three parameters and an opaque function.

# modelterms7 0.45.0

* A smooth with a factor `by` says its penalty is one copy per level rather
  than handing `penalties7` the assembled product: `quadratic_penalty(P,
  blocks = m)` instead of `quadratic_penalty(kronecker(diag(m), P))`. Nothing
  of size `(mk)^2` is decomposed, and at `m = 200` over a basis of ten the
  term builds in 0.06 s where it took 4.16 s.

* The ANISOTROPIC tensor branch still assembles: `additive_penalty()` reads
  one eigendecomposition of the SUM of its components, which is not a
  blockwise quantity, so the same shortcut does not apply to it.

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
  modeling statement. Nor could the value be policed in general: where the
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
