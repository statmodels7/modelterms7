# Check a Term's Grid Sizes Against Its Penalty

Validates the number of values that the constructor was given per
hyperparameter, and returns them as a named list from which the `NULL`
entries are dropped.

## Usage

``` r
check_grid(values, penalty, what = "this term")
```

## Arguments

- values:

  A named list of the constructor's arguments, `NULL` where the path's
  default size is wanted.

- penalty:

  A penalties7 penalty, or a function returning one, used only to read
  the names.

- what:

  The constructor's name, for the message.

## Value

A named list of grid sizes.

## Details

How finely a hyperparameter is swept belongs to the term for the same
reason as whether it is swept at all: a penalized block of four columns
and one of four hundred need different grids, and a criterion is applied
to every term of the model at once. Where a term gives no size, the path
of the fitting layer uses its own default.

The smallest grid has two values. There is no upper limit; each point of
a path is a whole fit, so the cost grows linearly in the number of
values on one axis.

## See also

[`check_hyper()`](https://statmodels7.github.io/modelterms7/reference/check_hyper.md),
[`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md)
