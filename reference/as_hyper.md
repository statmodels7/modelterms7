# Read a Caller's Held Hyperparameters

Accepts a named vector or a named list and returns a named list,
checking only what can be checked before the penalty exists.

## Usage

``` r
as_hyper(x, what = "this term")
```

## Arguments

- x:

  A named vector, a named list, or `NULL`.

- what:

  The term's label, for the message.

## Value

A named list, possibly empty.

## Details

[`random`](https://statmodels7.github.io/modelterms7/reference/random.md)
builds one of three penalties depending on what it was given – a ridge,
a structured prior over a parameters7 matrix, or a distributions7 family
used coordinatewise – so which names there are is not known until the
term is built. The shape is checked here and the names against the
penalty by
[`check_hyper`](https://statmodels7.github.io/modelterms7/reference/check_hyper.md)
at that point.

A vector entry is a written-out grid rather than a held value, exactly
as it is in the constructors that name their hyperparameters, so the
length is not checked here either;
[`.hyper_parts`](https://statmodels7.github.io/modelterms7/reference/dot-hyper_parts.md)
splits the two once the penalty exists.
