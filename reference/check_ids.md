# Check the Sharing Labels That a Constructor Was Given

Validates `id` and returns it as a named character vector, one entry per
hyperparameter to share, empty where the constructor was given none.

## Usage

``` r
check_ids(x, names, what = "this term")
```

## Arguments

- x:

  The `id` argument as given: `NULL`, a single string, or a named
  character vector.

- names:

  The hyperparameter names that the term's penalty carries, against
  which the names of the labels are checked, or `NULL` where the penalty
  does not exist yet. In that case the shape is checked and the labels
  are carried as written, to be resolved by a second call where the
  names are known.

- what:

  How to name the term in a message.

## Value

A character vector, possibly empty, named by the hyperparameters. In the
deferred case (`names = NULL`) a single unnamed string is returned
unnamed, to be named by the second call.

## Details

Two forms are accepted and both have the same meaning. An unnamed single
string shares the hyperparameter, and is valid only where the penalty
carries one: it is the short form for the frequent case, a smooth or a
ridge. A named vector, `c(alpha = "A")`, names the hyperparameter to
share, and a penalty carrying several requires this form, since sharing
the mixing of an elastic net while leaving its strength free is a
different model from the reverse, and the code does not infer which one
was meant.

The short form on a penalty carrying several hyperparameters is an error
that lists them. The names are checked against the penalty's own, so a
misspelling is reported where it was written.

Sharing does not merge the penalties. They remain two objects estimated
at one value, and the code does not check whether they are comparable:
the two may weigh quite different matrices, and whether that is the
intended model is left to the caller. A shared value is most meaningful
between penalties of the same shape over covariates of comparable scale.

## See also

[`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md),
which reports it;
[`check_hyper()`](https://statmodels7.github.io/modelterms7/reference/check_hyper.md)
for the argument that holds a hyperparameter instead of sharing it.
