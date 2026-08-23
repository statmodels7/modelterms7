# Continuing a Score-Driven Filter Past the Series

The level the recursion reaches at rows that come after the observed
ones.

## Arguments

- term:

  A built score-driven term.

- psi:

  The term's parameters, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- f_past:

  The level at each observed row.

- s_past:

  The score at each observed row.

- newdata:

  The rows to continue onto.

- ...:

  Ignored.

## Value

A numeric vector of `nrow(newdata)` levels.

## Details

The recursion is the filter's own, started from the level and score the
observed part ended at rather than from the term's own seed. What
changes past the data is the score: it has zero conditional mean by
construction, the model's own defining property, so at a row whose
response is not observed the driving term is its expectation and the
continuation is the DETERMINISTIC recursion \$\$f\_{n+h} = \omega +
\sum_i \alpha_i s\_{n+h-i} + \sum_j \beta_j f\_{n+h-j},\$\$ the loadings
contributing only while \\n+h-i\\ is still an observed time. No
simulation and no integration is involved.

A new row is placed by its own time within its own group, and must come
after every observed time of that group: a row falling inside the series
is not a continuation but a re-reading of it, where the response is
known and the filter must be run rather than continued, so it is
rejected with the rows named. A group the fit never saw is rejected for
the same reason – there is no state to continue.

## See also

[`term_continue`](https://statmodels7.github.io/modelterms7/reference/term_continue.md),
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
