# What a Grouping Variable Has to Be

Rejects a grouping expression whose value is continuous, naming it and
the spelling that says the values really are labels.

## Usage

``` r
.random_check_group(v, expr)
```

## Arguments

- v:

  The value of the grouping expression.

- expr:

  The expression, for the message.

## Value

`v`, unchanged; called for its error.

## Details

A grouping variable's values are labels, and
[`factor()`](https://rdrr.io/r/base/factor.html) turns any vector into
one by formatting its values, so a continuous covariate in that position
gives one level per distinct value and a random effect with as many
coefficients as there are rows. Measured before the check existed,
`random(~ 1 | x)` with a standard normal `x` at 60 observations built 60
columns and reported nothing.

What is rejected is a **double carrying a value that is not a whole
number**, which is narrow on purpose. Integer codes are the ordinary way
to write a grouping variable and stay legal, and a caller who means the
levels of a non-integer vector writes
[`factor()`](https://rdrr.io/r/base/factor.html) around it, which
arrives here as a factor and passes.

The check also does the work the two-bar form needs: `~ 1 | u | id` and
`~ 1 | id | u` are syntactically the same shape and differ only in which
position holds what, so the defence is that the last position must be
usable as a grouping variable.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
whose `formula` argument documents the rule;
[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md),
which reaches it through the grouping variable's evaluation.
