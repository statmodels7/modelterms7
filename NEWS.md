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
