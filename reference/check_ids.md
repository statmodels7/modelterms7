# Check the Sharing Labels a Constructor Was Given

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

  The hyperparameter names the term's penalty carries, which the labels'
  names are checked against, or `NULL` where the penalty does not exist
  yet: the shape is then checked and the labels carried as written, to
  be resolved by a second call where the names are known.

- what:

  How to name the term in a message.

## Value

A named character vector, possibly empty.

## Details

Two forms are accepted and they say the same thing. An unnamed single
string shares **the** hyperparameter, and is legal only where the
penalty carries one: it is the short form for the frequent case, a
smooth or a ridge. A named vector names which, `c(alpha = "A")`, and is
what a penalty carrying several needs, since sharing the mixing of an
elastic net while leaving its strength free is a model and guessing
which was meant is not.

The short form on a penalty carrying several is an error listing them,
rather than a rule about position. The names are checked against the
penalty's own, so a misspelling is reported where it was written.

Sharing does not merge the penalties. They stay two objects estimated at
one value, so nothing here asks whether they are comparable: the two may
weigh quite different matrices, and whether that is the model wanted is
the caller's to know. What the documentation says instead is that `id`
is worth care, and that it means what one expects between penalties of
the same shape over covariates of comparable scale.

## See also

[`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md),
which reports it;
[`check_hyper()`](https://statmodels7.github.io/modelterms7/reference/check_hyper.md)
for the argument that holds a hyperparameter instead of sharing it.
