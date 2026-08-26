# Generating From a Score-Driven Filter

Runs the recursion forward drawing the response at each step.

## Arguments

- term:

  A built score-driven term.

- psi:

  The term's parameters, on the parameter scale.

- eta:

  The static part of the predictor.

- draw:

  A function `(e, i)` returning one response value.

- score, curvature:

  Functions of `(y, e, i)` giving the derivatives of the log-density in
  the predictor at observation `i`. They take the response as an
  argument, unlike the filter's own, because here it is being made.

- ...:

  Ignored.

## Value

A list with `eta`, `y` and `latent`, the level.

## Details

No separate recursion is written. The filter's own recursion is the
generative one. The level at one time is a function of the scores before
it, and a score is a function of a response and a predictor, so the only
difference is where the response comes from.
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
calls its `score` callback exactly once per observation, in time order
within each group, at the predictor the recursion has just produced; a
callback that draws the response there, keeps it, and returns the score
of what it drew turns the filter into a generator.

The curvature is read at the same point, so it reads the response the
score drew, drawing no second one.

The fast route is not taken: it reads the response through a registered
C entry point, and here the response does not exist until the step that
needs it.

## See also

[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md),
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
