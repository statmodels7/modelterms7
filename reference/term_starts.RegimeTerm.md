# The Starting Points of a Regime Term

`n_start` starting points. The first is
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)'s.
Each further one draws every additive log-ratio of the transition matrix
from \\N(0, 2^2)\\ and adds \\N(0, 0.7^2)\\ to the logarithm of every
gap, the first level staying where the first start put it.

## Arguments

- term:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- ...:

  Passed to
  [`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md).

- target:

  The response on the scale of the predictor, or `NULL`.

## Value

A list of `term@n_start` named numeric vectors on the unconstrained
scale.

## Details

A log-ratio of standard deviation 2 reaches a transition probability of
0.02 or 0.98 about one time in three, so the draws cover chains that
stay in a regime and chains that leave it at once. A gap scaled by
\\e^{\pm 0.7}\\, about a factor of two either way, keeps the levels
ordered and within the range of the response. The draw for start \\s\\
is made with the seed \\100 + s - 1\\, and the caller's random number
generator is restored afterwards, so the same fit gives the same starts
and does not move the caller's stream. Measured on
[`MASS::geyser`](https://rdrr.io/pkg/MASS/man/geyser.html) with three
regimes, eight starts reach five distinct maxima of the log-likelihood,
from -1053.39 to -1210.49.
