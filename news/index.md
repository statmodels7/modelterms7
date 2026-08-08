# Changelog

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
