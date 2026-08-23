# Drawing a Regime Path

A path of the latent chain, and the predictor each observation gets from
the level of the state it lands in.

## Arguments

- term:

  A built regime term.

- psi:

  The term's parameters, on the parameter scale.

- eta:

  The static part of the predictor.

- draw:

  Ignored: the levels do not read the response.

- ...:

  Ignored.

## Value

A list with `eta`, `y` (`NULL`) and `latent`, the state of each
observation.

## Details

The chain is started from its STATIONARY distribution rather than from
the first state or from a uniform draw: the model's likelihood is
written with that initial law, so any other start would simulate a
different model from the one a fit reads back. The transition matrix and
the stationary law come from the term's own
[`transition_matrix`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html),
so nothing about the chart is restated here.

The response is not drawn – the levels do not read it – so the caller
draws at the returned predictor, and `y` is `NULL`.

## See also

[`term_simulate`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md),
[`regime`](https://statmodels7.github.io/modelterms7/reference/regime.md)
