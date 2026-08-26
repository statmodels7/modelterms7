# Check a Term's Held Hyperparameters Against Its Penalty

Validates the values a constructor was given against the names and the
bounds of the penalty the term will build, and returns them as a named
list, dropping the `NULL` entries, which are the ones to be estimated.

## Usage

``` r
check_hyper(values, penalty, what = "this term")
```

## Arguments

- values:

  A named list of the constructor's arguments, `NULL` where the
  hyperparameter is to be estimated.

- penalty:

  A penalties7 penalty, or a function returning one, used only to read
  the names and the bounds.

- what:

  The constructor's name, for the message.

## Value

A named list, keyed by the penalty's own hyperparameter names, of the
values given exactly one number. Empty where every hyperparameter was
left `NULL` or written out as a grid.

## Details

The check happens at construction, where the caller can see it, instead
of at the fit three layers away. A name the penalty does not carry is an
error naming the ones it does. That is what turns `mcp(x, a = 3)`,
SCAD's shape written on an MCP whose own shape is `gamma`, into a
message instead of an argument that lands in `...` and does nothing.

The bounds are open, as they are everywhere in the toolkit: a ridge at
\\\lambda = 0\\ is no penalty at all, and an elastic net at \\\alpha =
0\\ has no kink and is a penalty of another kind.

One argument carries three states, each read for one hyperparameter at a
time: `NULL` has the path build the grid, one number holds the
hyperparameter, and several are the grid itself. This returns the second
of the three;
[`check_values()`](https://statmodels7.github.io/modelterms7/reference/check_values.md)
returns the third.

## See also

[`check_values()`](https://statmodels7.github.io/modelterms7/reference/check_values.md),
[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md),
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
