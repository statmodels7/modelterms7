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
  argument, unlike the filter's own, because the response is generated
  during the pass. The defaults, used where the argument is `NULL`, are
  `y - e` and `-1` (a Gaussian response with unit variance).

- ...:

  Ignored.

## Value

A list with `eta`, `y` and `latent`, the level.

## Details

The generative recursion is the filter's own. The level at one time is a
function of the scores before it, and a score is a function of a
response and a predictor, so the filter and the generator differ only in
where the response comes from.
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
calls its `score` callback exactly once per observation, in time order
within each group, at the predictor that the recursion has just
produced. A callback that draws the response at that point, stores it
and returns the score of the drawn value turns the filter into a
generator.

The curvature is read at the same point, so it uses the response that
the score callback has already drawn and draws no second one.

The fast route is not taken: it reads the response through a registered
C entry point, and here the response does not exist until the step that
needs it.

## See also

[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md),
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
