# Check a Term's Written-Out Grids Against Its Penalty

Validates the hyperparameter values a constructor was given as a vector,
which a path visits as they stand.

## Usage

``` r
check_values(values, penalty, what = "this term")
```

## Arguments

- values:

  A named list of the constructor's arguments, `NULL` where the
  hyperparameter is to be estimated.

- penalty:

  A penalties7 penalty, or a function returning one, used only to read
  the names and the bounds.

- what:

  The constructor's name, for the message.

## Value

A named list of numeric vectors, one entry per hyperparameter the caller
wrote out.

## Details

Several values are neither a held hyperparameter nor a request to build
a grid: they are the grid. Each is checked against the penalty's bounds
as a held value would be. Nothing else is applied to them: the value
that empties the block does not cap them and the depth of the path does
not reach them, the caller having said which values to visit.

They are sorted, because a path is walked from the emptiest fit towards
the fullest and its warm starts follow that order, and duplicates are
dropped. Which end of the order is the sparse one depends on the
penalty, so the direction is settled where the path is built.

## See also

[`check_hyper()`](https://statmodels7.github.io/modelterms7/reference/check_hyper.md),
[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
