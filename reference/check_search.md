# Check How a Term Covers Its Own Hyperparameters

Validates the sweep a constructor was given for the term's own kinked
hyperparameters.

## Usage

``` r
check_search(v, what = "this term")
```

## Arguments

- v:

  What the constructor was given, or `NULL`.

- what:

  The term's label, for the message.

## Value

A character vector of length one, or of length zero.

## Details

The sweep is one word per term, because it states how the
hyperparameters are combined with each other, which is not a property of
any single hyperparameter. It belongs to the term because a penalty with
a kink is fitted by a scheme of its own, and the way that scheme sweeps
its hyperparameters is part of it. A criterion is applied to every
hyperparameter of the model, the smooth ones included, and most of them
cannot use such an argument.

## See also

[`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md),
[`check_min_ratio()`](https://statmodels7.github.io/modelterms7/reference/check_min_ratio.md)
