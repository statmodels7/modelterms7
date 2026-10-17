# Drawing Break-Points From Their Prior

Draws one set of latent positions per group from the prior that the term
declares, and returns the predictor that each observation receives from
them.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- psi:

  The term's parameters, on the parameter scale.

- eta:

  The static part of the predictor.

- draw:

  Ignored: the positions do not depend on the response.

- ...:

  Ignored.

## Value

A list with `eta`, `y` (`NULL`) and `latent`, a data frame of the drawn
positions with the columns `group`, `psi` (the break-point) and `value`.

## Details

The latent positions are integrated out of the likelihood and never
estimated, so simulating from the model consists of drawing them, once
per group, and then evaluating the term at the drawn positions. Under
the gaussian prior the draw is \\N(m_k, \tau_k^2)\\, where \\m_k\\ is
`psik.mean` and \\\tau_k\\ is `psik.sd`. Under an explicit prior it is a
draw from that family, whose location is fixed at zero, shifted by
`psi1.mean`, which is the convention of the likelihood.

The shift is the term's own construction evaluated at the drawn
positions: a change of level at each break-point for the step kind, a
change of slope for the continuous kind, both for the joint kind, and
the linear term beside them where the term carries it.

The response is not drawn, because the positions do not depend on it, so
the caller draws the response at the returned predictor.

## See also

[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
