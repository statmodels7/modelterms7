# The Grouping of a Random-Effect Term

The grouping expression, its levels and the number of within-group
columns, or `NULL` on an unbuilt term.

## Arguments

- term:

  A built
  [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md).

- ...:

  Unused.

## Value

A list with `expr`, `levels` and `dim`, or `NULL`.

## Details

The three come from the blueprint, where the build recorded them. `dim`
is the within-group column count, so the block's level \\i\\ occupies
positions \\(i-1)d + 1\\ to \\id\\, which is the group-major order the
blockwise penalty reads.

## See also

[`term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.md)
for the generic.
