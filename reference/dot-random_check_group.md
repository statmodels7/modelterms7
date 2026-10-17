# What a Grouping Variable Has to Be

Rejects a grouping expression whose value is continuous, naming it and
the spelling to use when the values really are labels.

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
would give one level per distinct value and a random effect with as many
coefficients as there are rows.

Only a double carrying a value that is not a whole number is rejected.
Integer codes are the ordinary way to write a grouping variable and are
accepted, and a caller who means the levels of a non-integer vector
writes [`factor()`](https://rdrr.io/r/base/factor.html) around it, which
arrives here as a factor and passes.

The check also applies to the two-bar form: `~ 1 | u | id` and
`~ 1 | id | u` have the same syntactic shape and differ only in which
position holds which element, so the last position is required to be
usable as a grouping variable.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
whose `formula` argument documents the rule;
[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md),
which reaches it through the grouping variable's evaluation.
