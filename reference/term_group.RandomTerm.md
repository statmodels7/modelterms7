# The Grouping of a Random-Effect Term

Returns the grouping expression, its levels, the number of within-group
columns and their names, or `NULL` on an unbuilt term.

## Arguments

- term:

  A built
  [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md).

- ...:

  Unused.

## Value

A list with `expr`, `levels`, `dim` and `names`, or `NULL`.

## Details

These come from the `build_info`, where the build recorded them. `dim`
is the within-group column count, so the block's level \\i\\ occupies
positions \\(i-1)d + 1\\ to \\id\\, which is the group-major order that
the blockwise penalty reads. `names` holds the column names of the
within-group design.

## See also

[`term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.md)
for the generic.
