# Generating From a Joint Score-Driven Filter

Runs the joint recursion forward, drawing the response at each step.

## Arguments

- term:

  A built
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- psi:

  The parameters, as for
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- eta:

  The static part of the filtered predictors, a matrix with one column
  per filtered parameter.

- draw:

  A function `(e, i)` of the vector of filtered predictors at row `i`,
  returning one response value.

- score, curvature:

  Functions `(y, e, i)` returning the driving vector and its derivative
  in the filtered predictors, as the callbacks of
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
  do, at the response `y`.

- ...:

  Ignored.

## Value

A list with `eta` (the filtered predictors), `y` and `latent` (the
levels), matrices with one column per filtered parameter apart from `y`.

## Details

As for one term
([`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md)):
the filter calls its driving vector once per observation, in time order
within each group, at the filtered predictors that it has just produced,
so a callback that draws the response there, stores it and returns the
driving vector of the drawn value turns the filter into a generator. The
curvature is read at the response already drawn.

## See also

[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md),
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
