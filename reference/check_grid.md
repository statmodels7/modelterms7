# Check a Term's Grid Sizes Against Its Penalty

Validates the number of values the constructor was given per
hyperparameter, and returns them as a named list with the `NULL` entries
dropped.

## Usage

``` r
check_grid(values, penalty, what = "this term")
```

## Arguments

- values:

  A named list of the constructor's arguments, `NULL` where the
  criterion's default is wanted.

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
and one of four hundred want different grids, and a criterion applies to
every term of the model at once and cannot know which it is looking at.
Where a term says nothing the criterion's own default is used.

Two values is the smallest grid that is a grid. There is no upper limit
beyond the caller's patience: each point of a path is a whole fit, so
the cost is linear in the number asked for.

## See also

[`check_hyper()`](https://statmodels7.github.io/modelterms7/reference/check_hyper.md),
[`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md)
