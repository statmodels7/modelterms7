# The Filter's Parameters at Rows Outside the Fitting Data

A developed parameter's value at each new row, read through each
sub-term's own blueprint rather than rebuilt.

## Usage

``` r
.gas_sub_new(term, u, newdata)
```

## Arguments

- term:

  A built score-driven term.

- u:

  The term's free vector.

- newdata:

  The rows to read at.

## Value

A list with `om`, `A` and `B`, one row each per observation of
`newdata`.

## Details

It is
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
on every sub-term, which is what keeps a basis or a set of contrasts
from being relearned at other rows, and then the same chart the fit
used. The persistence goes through Levinson-Durbin exactly as it does
inside the filter.

## See also

[`term_continue`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)
