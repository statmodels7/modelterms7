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

One word per term, one per hyperparameter being meaningless: it says how
the hyperparameters are combined with each other, which is not a
property any one of them has. It belongs to the term because a penalty
with a kink is fitted by a scheme of its own, and how that scheme sweeps
its own hyperparameters is part of the scheme. A criterion applies to
every hyperparameter of the model, the smooth ones included, and would
be carrying an argument most of them cannot read.

## See also

[`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md),
[`check_min_ratio()`](https://statmodels7.github.io/modelterms7/reference/check_min_ratio.md)
