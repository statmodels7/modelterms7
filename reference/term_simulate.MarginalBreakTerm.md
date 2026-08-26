# Drawing Break-Points From Their Prior

One set of latent positions per group, drawn from the prior the term
declares, and the predictor each observation gets from them.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- psi:

  The term's parameters, on the parameter scale.

- eta:

  The static part of the predictor.

- draw:

  Ignored: the positions do not read the response.

- ...:

  Ignored.

## Value

A list with `eta`, `y` (`NULL`) and `latent`, a data frame of the drawn
positions by group.

## Details

The latent positions **are** the model here, integrated out of the
likelihood and never estimated, so simulating from the model means
drawing them, once per group, and then evaluating the term at what was
drawn. Under the gaussian prior that is \\N(m_k, \tau_k)\\; under an
explicit prior it is a draw from that family with its location fixed at
zero, shifted by \\m_1\\, which is the same convention the likelihood is
written with.

The shift is the term's own construction read at the drawn positions: a
change of level at each break-point for the step kind, a change of slope
for the continuous one, both for the joint one, and the linear term
beside them where the term carries it.

The response is not drawn, the positions not reading it, so the caller
draws at the returned predictor.

## See also

[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
