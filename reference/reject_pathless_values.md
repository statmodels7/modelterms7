# Reject a Written-Out Grid Where There Is No Path

Signals an error when the caller gave several values for a
hyperparameter whose penalty has no kink.

## Usage

``` r
reject_pathless_values(values, pen, what = "this term")
```

## Arguments

- values:

  The written-out grids, as
  [`check_values()`](https://statmodels7.github.io/modelterms7/reference/check_values.md)
  returns them.

- pen:

  The penalty the term will build, or `NULL`.

- what:

  The term's label, for the message.

## Value

`NULL`, invisibly.

## Details

Several values are a grid for a path to visit, and only a penalty with a
kink is swept along one: everything else has its hyperparameter read at
the mode by a marginal criterion, which would take the vector and do
nothing with it. The question is put to the penalty at a probe value of
its own hyperparameters, the kink set being structural, so a ridge and a
random effect under a Gaussian prior are covered by the same line as a
family added later.

## See also

[`check_values()`](https://statmodels7.github.io/modelterms7/reference/check_values.md),
[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
