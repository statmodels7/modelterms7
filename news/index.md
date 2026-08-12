# Changelog

## modelterms7 0.20.0

- [`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
  returns the exact Hessian of a likelihood mixed over latent states, in
  the whole of a caller’s unknown vector: the coefficients of every
  equation together with the term’s own parameters.

  What a caller could assemble from
  [`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
  alone is the COMPLETE-DATA information, the ordinary one averaged over
  the smoothed states. That is the matrix an EM step inverts and it is
  not the observed information: by the missing-information principle the
  two differ by the conditional variance of the complete-data score, so
  the complete-data one is the larger and a standard error read off it
  is too small. Measured on a two-regime gaussian, the difference
  reaches 30 per cent of a standard error where the regimes overlap and
  vanishes as they separate and the states become known.

  Louis’s identity is one route to it and is not the one taken. The
  scaled forward recursion computes the observed log-likelihood exactly,
  as a sum of the logarithms of its normalizing constants, so
  differentiating that arithmetic twice gives the observed Hessian with
  no identity, no pairwise smoothed probabilities and no second-moment
  recursion. The first and second derivatives of the filtered
  distribution are propagated beside it and renormalized by the quotient
  rule. Louis’s identity becomes the check instead: the difference
  between the two matrices must be positive semidefinite, and it is
  measured to be, strictly, wherever the states carry any uncertainty.

  The cost is `O(n K^2 m^2)`, and the computation is meant to run once
  at a fitted point.
  [`regime_stationary()`](https://statmodels7.github.io/modelterms7/reference/regime_stationary.md)
  gained the second derivative of the stationary distribution, from the
  same linear system as its first.

## modelterms7 0.19.0

- [`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
  returns the smoothed state probabilities of a latent Markov term,
  which is everything a model layer needs to differentiate a likelihood
  mixed over states. By Fisher’s identity the derivative of that
  likelihood in ANY predictor the model carries is the
  posterior-weighted derivative of the ordinary one, so a caller
  differentiates its own log-density K times vectorized and weights, and
  needs no callback per observation. That is the property that made the
  forward pass compilable read once more: a regime shifts a predictor
  known before the recursion starts.

  The probabilities come from the forward pass this term already runs
  and a backward pass beside it, both normalized – without which the
  quantities are products of t densities and reach zero in double
  precision within a few hundred observations. Validated against
  `numDeriv`: the rows sum to one to 1.1e-16 and Fisher’s identity holds
  to 8.1e-09, 1.3e-07 and 7.7e-10 over one series, three regimes and
  groups, where the score at the marginal mean is out by 1.4.

- [`term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)
  says which of a term’s parameters shifts its equation’s predictor by a
  constant: `"omega"` for
  [`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
  `"level1"` for
  [`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md),
  and `character(0)` for everything else. It exists so that a fitting
  layer can resolve the confounding with an intercept rather than refuse
  the model. Which parameter is the level is the term’s answer; which
  one is dropped is the layer’s, since only the layer knows what else
  the equation carries.

## modelterms7 0.18.0

- [`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
  is the second-order companion of
  [`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md):
  the forward Jacobian of the predictor a structural term produces in a
  caller’s unknowns, and the second derivative of that predictor
  contracted against the caller’s weights. It is what an observed
  information needs and what the reverse recursion alone does not give.

  The contract keeps the split the adjoint already uses. `seed` is the
  derivative of the static predictor in the caller’s unknowns, so the
  term learns nothing else about them; `blocks` is a callback returning,
  at an observation and the Jacobian the recursion has reached, the two
  model quantities that seed the first and second derivatives of the
  score – `sum_q l_pq C_q` and `sum_{r,r'} l_prr' V_r' V_r'`. A model of
  one equation supplies zero and `l_ppp D'D`.

  Against numDeriv at p, q in {1,2}^2: the Jacobian 1.05e-10 and the
  contracted second derivative 1.0e-09 on a matrix of scale 24.6.

  Deviations are refused rather than silently mishandled: the per-group
  chain adds a factor to every derivative and is not written.

- [`gas_levinson2()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson2.md)
  carries the second derivatives of the autoregressive coefficients in
  the partial autocorrelations, which the curvature needs because the
  persistence reaches the predictor through that map. The recursion is
  bilinear, so differentiating twice adds no new kind of term, only the
  two places the product rule puts the first derivative. Against
  numDeriv: 3.8e-12, 2.9e-12 and 9.4e-11 at q = 2, 3 and 4. The value
  and the jacobian are bit-identical to
  [`gas_levinson()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson.md),
  which is what says this is not a second route to them, and the last
  coefficient’s second derivative is exactly zero at every order, it
  being the last partial autocorrelation.

## modelterms7 0.17.0

- [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)’s
  numerical route uses numericals7’s stencil library. It wrote out a
  three-point central difference with its own step, and modelterms7 did
  not import numericals7 at all – the nodes, the weights and the step
  are `fd_offsets()`, `fd_weights()` and `fd_step()` now, at accuracy
  four, which is the five-point rule. Measured against the exact
  symbolic Jacobian of the same function, on the route an opaque
  `f(x, theta)` takes: 3.7e-13 against the three-point rule’s 2.7e-11.
  This Jacobian is the design block, so its accuracy is the accuracy of
  every step such a fit takes; the two extra evaluations are the trade
  distributions7 measured for the skew t.

- [`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
  is the reverse recursion of a structural term: the derivative of a
  caller’s objective with respect to the static predictor the term was
  handed, and with respect to the sequence of scores it was given.

  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
  returns the derivative of the predictor in the term’s OWN parameters,
  which is what estimating those needs, and it is not what estimating
  the coefficients of the same equation needs. A score-driven level at
  one time is driven by the scores at earlier ones, read at predictors
  those coefficients also enter, so the derivative of the predictor in a
  coefficient carries a term the block does not. Measured against
  `numDeriv` on the derivative in the static predictor: the reverse
  recursion agrees to 1e-8 and the direct score alone is wrong by 0.6 to
  1.05 in every configuration tried – one series, p = 2 and q = 2,
  groups, and groups with deviations.

  Propagating that forward would cost one derivative array per
  coefficient; the reverse pass costs one whatever their number. Two
  derivatives are returned rather than one because the score depends on
  more than the predictor it is read at: multiplying `dscore` by the
  mixed second derivative of a log-density gives the derivative in
  ANOTHER equation’s predictor, which is what a model layer with several
  distribution parameters needs.

- [`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
  takes `newdata`, so a fitting layer can compute a term’s contribution
  on other rows. Where the block is a Jacobian,
  [`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  times the coefficients is the linearization and not the contribution:
  for
  [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
  the two differ by a step at the break-point in a construction that is
  continuous. Rows are treated as
  [`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  treats them, through the levels and constants the blueprint recorded.

- [`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
  asks whether a term’s own iteration has settled, which a score cannot
  always answer. Where the block is the Jacobian of the contribution the
  gradient of the model’s objective is the model’s and its vanishing is
  the test; where the block is a working linearization with a frozen
  weight, as in
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
  and
  [`jseg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
  the profile objective is a step function in the break-point and has no
  gradient to vanish. The base method is `TRUE`, and the segmented
  method is
  [`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md).

## modelterms7 0.16.0

- [`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
  carries a population value and a deviation per group. With `by` and
  `deviations` each group of a panel is filtered with parameters of its
  own, written as

  ``` R
  psi[j, i] = g_j^-1( g_j(psi_j) + delta[j, i] ),
  ```

  the deviation acting on the unconstrained scale of the chart the
  parameter lives on, so a persistence stays inside (-1, 1) whatever the
  deviation is. `deviations` takes TRUE for every parameter or the names
  of the ones that carry one, and needs `by`. The deviations are
  parameters of the term, named after the parameter and the level
  (`omega.dev.a`), and carry the identity link, being unconstrained
  already.

  They are parameters and NOT a penalty on the per-group values through
  a difference matrix, which is the same model written the other way.
  The difference decides what can be fitted: a penalty over a general
  map is the generalized-lasso problem, whose proximal operator does not
  split by coordinate, while a deviation named as a coordinate is
  reached by a soft threshold and by a coordinate descent unchanged.

  The filter runs once per group and chains the columns of its jacobian
  onto the population values, exactly, `d psi[j,i] / d psi_j` being
  `g^-1'(g(psi_j) + delta) g'(psi_j)` and the derivative in the
  deviation the same without the second factor. At a zero deviation the
  two are reciprocal by the inverse function theorem, so the filter and
  every population column are then bit for bit the shared-parameter
  ones. The jacobian agrees with `numDeriv` to 1e-10 with deviations on
  one parameter and on all of them.

  The deviations are identified by their penalty and not otherwise. A
  parameter and its m deviations are m+1 numbers describing m group
  values, so a constant added to the population value on the
  unconstrained scale and subtracted from every deviation leaves the
  filter exactly unchanged, and the likelihood is flat along one
  direction per parameter carrying them. That is the parametrization of
  a random effect, identified there by a variance component and here by
  the penalty, which selects the deviations of smallest size among the
  descriptions of the same panel.

- `gas(penalty =)`, `nl(penalty =, penalize =)` and `seg(penalty =)`
  declare the parameters they penalize through
  [`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
  naming the coordinates each penalty covers instead of selecting them
  from the block with a map. The map was the defect: a separable penalty
  under a selection map is the generalized lasso, so
  [`penalties7::has_prox()`](https://statmodels7.github.io/penalties7/reference/has_prox.html)
  was FALSE for a `seg(penalty = "lasso")` and neither a proximal step
  nor a coordinate descent could be taken on it. Named as coordinates
  the map is the identity and both are available.

  - [`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
    penalizes the deviations, one penalty per parameter carrying them,
    and rejects a penalty without them: the population parameters of a
    filter are not shrunk towards zero.
  - [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
    penalizes one parameter at a time, the whole coefficient vector
    where the parameter carries a subformula, so a lasso there selects
    which covariates a parameter depends on. What is shrunk is the
    coefficient, so with a link the target is `g^-1(0)` and not zero.
  - [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
    penalizes the changes as before, and
    [`jseg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
    now declares two penalties rather than one over their union: a slope
    change and a jump are not comparable quantities and cannot share a
    hyperparameter.

- [`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
  counts a term parameter by parameter rather than reading one penalty
  for the whole block. A parameter no penalty reaches counts one; a
  parameter under a kinked penalty counts one when it is away from zero;
  the rest are counted together by `tr[(H+S)^-1 H]` over the sub-block
  they occupy, with `S` carrying each smooth penalty’s Hessian at the
  parameters it covers and zero elsewhere. Each rule reduces to what the
  term reported before when one penalty covers the whole block. `theta`
  is that penalty’s hyperparameters for a term carrying one, and a list
  keyed by the penalty names for a term carrying several.

- [`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
  reads
  [`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
  too, so a term penalized over part of its parameters answers for the
  part: `seg(x, penalty = "lasso")` is not smooth although its linear
  effect and its break-points are unpenalized.

- [`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
  answers for a structural term as well, counting the entries of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
  which is the vector
  [`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
  indexes into there.

## modelterms7 0.15.0

- [`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
  is what a term declares it wants penalized: a list of entries, each
  naming a subset of the term’s own parameters and the penalty over
  them. It generalizes
  [`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
  in the two directions a model layer needs.

  A term may carry MORE THAN ONE penalty, over different parameters of
  its own – a panel model with a population value free and a deviation
  per group shrunk is one penalty over part of the parameters and none
  over the rest. And the parameters need NOT be coefficients of a design
  block: the persistence of a score-driven term, the nonlinear
  parameters of
  [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
  the break-point of
  [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
  are parameters of the term and nothing else, and all a penalty needs
  of them is a vector of numbers and their positions.

  The base method answers from
  [`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
  so every term shipped here needs no method of its own and behaves
  exactly as before, and a structural term answers with an empty list
  rather than raising, which is what lets a caller enumerate over every
  term without knowing which kind it has.

  The entry’s name is unique WITHIN the term and is empty for a penalty
  over the whole of it. It is not the term’s name: two
  [`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  terms in one formula are two terms with their own hyperparameters, and
  it is the caller that knows what it called each one.

## modelterms7 0.14.0

- te() centers its block: the tensor product of the marginal bases
  contains the constant, and the null space of the tensor penalty
  contains it too, so beside an intercept the design was rank deficient
  by exactly one with no penalty covering the deficiency. Measured on
  300 observations of te(a, b, k = 5), the design of y ~ te(a, b, k = 5)
  had 25 of 26 columns, a smallest singular value of 2.2e-15 and a
  condition number of 8.0e15, and the penalized information a smallest
  eigenvalue at the rounding floor – which chol() accepts or rejects by
  the luck of rounding, so vcov(), confint() and the outer criterion
  returned numbers computed on a singular matrix. The block now carries
  the sum-to-zero constraint over the observed covariates, through
  basis7::constrain_basis() applied to the product basis, as mgcv does
  for a smooth. The same design is 25 of 25 columns at a condition
  number of 138.5, and the penalized information has a smallest
  eigenvalue of 0.93. A tensor term has one column fewer than the
  product of its marginal dimensions: te(x, z, k = 4) reports 15
  parameters where it reported
  16. The penalty follows by congruence and its rank does not move, the
      direction removed having been one of its null directions (21 of 24
      where it was 21 of 25). The transform is stored in the blueprint
      and reapplied by term_predict(), as the Demmler-Reinsch transform
      of s() is. The level of the surface is the model’s intercept, so a
      formula removing it fits a surface constrained to average zero.

## modelterms7 0.13.0

- enet(): the elastic-net term, beside ridge, lasso, scad and mcp,
  carrying penalties7::elasticnet_penalty(). Like the lasso it is not
  smooth, and its effective degrees of freedom are the nonzero count.

## modelterms7 0.12.0

- regime()’s forward recursion is compiled (src/regime_forward.cpp), and
  the density and score of every observation under every regime are
  computed once by k vectorized calls instead of 2nk scalar ones. Unlike
  the score-driven filter, nothing has to call back into R: a regime
  shifts the predictor by a level of its own, so none of those values
  depends on the filtered state. Measured against the R form, kept as
  the twin .regime_forward_r: 4.5x at k = 5, 13x at k = 3 and 28x at k =
  2, over T from 1e3 to 1e5, agreeing to 1.8e-15. End to end the term
  costs 2.04 microseconds per observation at k = 3, against 40 before.
- term_loglik()’s closures are called with the whole index vector. A
  closure returning one value where n were asked is rejected with a
  message saying so.

## modelterms7 0.11.1

- interpret_formula() rejects a call that evaluates to neither a model
  term nor a covariate, naming the call, its class and – when the
  function that was called is not the one modelterms7 exports under that
  name – the package that masked it. mgcv exports s() and te() and
  segmented exports seg(), so a user with either attached wrote our
  formula and got theirs; the value used to travel to model.matrix and
  fail there, naming neither the call nor the mask.

## modelterms7 0.11.0

- The discontinuous terms follow Fasola, Muggeo and Kuchenhoff’s
  algorithm as published. The weight is no longer capped: the covariate
  is rescaled away from the break-point by a factor `c0`, which leaves a
  gap around it, and
  [`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
  halves the factor whenever the break-point reverses direction. `band`
  is gone and no damping is needed, the factor governing the step.
  Measured over eight samples,
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
  now recovers the break-point from every starting position tried, where
  the capped form recovered it only from within a narrow basin.
- seg_step() and seg_converged() report the progress of the iteration
  and Fasola et al.’s stopping rule.
- seg_start() chooses the starting positions by scoring an equally
  spaced grid on the least-squares profile, which is the initialization
  Fasola et al. recommend and the piece that was missing. It is what
  settles the joint term: measured over eight samples of a jump and a
  change of slope at the same point, a single conventional start
  recovers the break-point in none to half of them depending on where it
  is placed, and the grid in all of them. Bootstrap restarting was
  measured beside it and does far less (0.12 to 0.75).
- Measured against the segmented package, on one covariate with one and
  two break-points and n from 200 to 20000: the continuous case agrees
  to four significant figures on the residual sum of squares and runs
  2.1 to 5.5 times faster in 4 to 7 iterations; the discontinuous case
  agrees where both start inside the basin and runs 1.1 to 1.8 times
  faster.

## modelterms7 0.10.0

- seg(), jump() and jseg(): break-points estimated with everything else.
  The continuous case is the Jacobian of Muggeo (2003); the
  discontinuous one uses the identity of Fasola, Muggeo and Kuchenhoff
  (2018), which is linear in the break-point once the weight is frozen,
  so the break-point is read off two coefficients rather than searched
  for. Both run under the existing term_refresh() contract, with by and
  an optional penalty on the changes. The working block is compiled
  (1.2x to 3.2x the R form over n from 1e3 to 1e6, agreeing to a
  rounding); the linear fit around it is BLAS in either language.

## modelterms7 0.9.0

- regime(): a latent Markov chain of regimes, each shifting the
  predictor by a level of its own, with the likelihood evaluated by the
  normalized forward recursion and its derivative propagated beside the
  state. Built on parameters7::transition_matrix(), which had no
  consumer until now.
- term_loglik(): the second shape of the structural branch, for a term
  whose contribution is a likelihood rather than a predictor.

## modelterms7 0.8.0

- nl(): the nonlinear parametric term. The design block is the Jacobian
  of the contribution in its parameters, refreshed by term_refresh() as
  they move, and term_value() reports the contribution a Gauss-Newton
  step needs beside it. The function may be a formula, read symbolically
  where deriv() manages it and differenced where it does not, or an
  opaque function, always differenced. Links per parameter; covariate
  submodels on the formula route, which is the only one that says where
  a parameter enters.

## modelterms7 0.7.0

- The score-driven recursion is compiled. The two callbacks into R
  remain, the score and the curvature belonging to the model’s
  distribution, but the arithmetic around them was 73 to 83 per cent of
  the loop’s time: measured 2.3x to 3.2x faster, and the R loop stays as
  the twin the kernel is compared against.

## modelterms7 0.6.0

- The structural branch is real: term_params(), term_links() and
  term_filter() define what a term that rewrites the likelihood must
  provide, and gas() implements score-driven dynamics over groups and
  time, with the persistence carried on a partial-autocorrelation chart
  and the derivative of the filter propagated alongside its state.
- s() and te(): penalized smooths of one and several covariates, the
  first under the Demmler-Reinsch reparametrization that separates the
  linear effect from the nonlinear deviation, both accepting `by`.

## modelterms7 0.5.0

- random() accepts slopes (~ x \| g) beside intercepts, with the
  within-group Gaussian unstructured or diagonal by `correlated`, a
  per-group precision structure replicated across groups through
  parameters7::kron_identity(), or a distribution applied
  coordinatewise.

## modelterms7 0.4.0

- random(~ 1 \| g): grouped random intercepts with the effect
  distribution as the penalty – independent Gaussian by default, a
  parameters7 precision structure, or a distributions7 object applied
  coordinatewise. Random slopes are rejected pending the block-diagonal
  composition in parameters7.

## modelterms7 0.3.0

- edf() with the counting rule per penalty (exact count, trace of the
  penalized smoother block, nonzero count), the penalty shown by print()
  on a built penalized term, and plot() at fitted coefficients.

## modelterms7 0.2.0

- The penalized terms ridge(), lasso(), scad() and mcp(), over formula
  or matrix input, with the penalty attached at build time and the
  smoothness flag read from its kink set.

## modelterms7 0.1.0

- First release: the term classes and generics, the unpenalized
  parametric term, the formula interpreter with recognition by
  evaluation, the censored-response constructor, and the term validator.
