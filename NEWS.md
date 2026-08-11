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
