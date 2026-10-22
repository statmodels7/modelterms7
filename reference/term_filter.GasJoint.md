# Filter Several Score-Driven Terms Jointly

Runs the joint recursion over each group in time order and returns the
filtered predictors together with their exact derivatives in the joint
parameter vector, propagated alongside the state.

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- eta:

  The static part of the predictors of the filtered parameters, a matrix
  with one column per entry of `term@params`.

- y:

  The response, unused directly.

- score:

  A function `function(e, i)` of the vector `e` of the filtered
  predictors at row `i`, returning the driving vector \\u_i\\, one value
  per filtered parameter.

- curvature:

  A function `function(e, i)` returning the matrix \\\partial
  u_i/\partial e^\top\\, whose entry \\(l, k)\\ is the derivative of the
  driving score of parameter \\l\\ in the predictor of parameter \\k\\.

- psi:

  A named list with one entry per filtered parameter, each the term's
  parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  of the term; or one flat list named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  of the joint filter.

- ...:

  Unused.

- fast:

  `NULL`, or the context of the compiled route, a list with `family`
  (the name under which the C registry of distributions7 knows the
  distribution), `constants` (its constants, each of length 1 or one per
  observation), `links` (the class name of each parameter's link),
  `link_par` (a list of each link's own parameters), `bounds` (a list of
  each link's bounds), `y`, `eta` (the static predictors of every
  distribution parameter, one column each), `filtered` (the positions of
  the filtered parameters among them) and `scaling`. Where the joint
  registry of distributions7 and the registry of linkfunctions7 cover
  the family and the links, the recursion runs compiled and reads the
  family through them instead of the callbacks; where they do not, the
  context is inert and the callbacks run in R. A scaled filter also
  needs the family's expected information and its derivative in the
  registry, which some families do not have (their expected information
  is a quadrature); their scaled filters run in R.

- threads:

  How many threads the compiled recursion may use, over groups, when the
  distribution's entries are safe on a worker thread and there are at
  least eight groups. Each group writes only to its own rows, so the
  result does not depend on the count.

## Value

A list with `eta` (a matrix like the input), `jacobian` (a named list
with one matrix per filtered parameter, one row per observation and one
column per entry of
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
of the joint filter) and `curv` (an array of dimension `c(n, K, K)`
holding the curvature read at each row).

## Details

With \\F\_{k,t}\\ the derivative of the level of parameter \\k\\ and
\\C_t = \partial u_t/\partial\eta\_{G,t}^\top\\ the derivative of the
driving vector, the derivative of the driving score of parameter \\l\\
is \\\sum_k C\_{t,lk} F\_{k,t}\\, because the static predictors do not
depend on the terms' parameters. The recursion of \\F\\ is the recursion
of the levels differentiated term by term.
