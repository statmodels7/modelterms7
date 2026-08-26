# Split a Constructor's Hyperparameter Arguments by Length

Validates each against the penalty and files it as held or as a grid.

## Usage

``` r
.hyper_parts(values, penalty, what = "this term")
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

A list of two named lists: `hyper`, the hyperparameters given one value
each, and `values`, those given several, sorted and deduplicated. Both
are keyed by the penalty's own hyperparameter names, and either may be
empty.

## Details

The validation is one body because the two states differ in length
alone: the name must be one the penalty carries and every value must lie
strictly inside its bounds, whether there is one of them or twenty.
