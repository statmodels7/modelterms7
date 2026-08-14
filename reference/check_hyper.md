# Check a Term's Held Hyperparameters Against Its Penalty

Validates the values a constructor was given against the names and the
bounds of the penalty the term will build, and returns them as a named
list with the `NULL` entries – the ones to be estimated – dropped.

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

A named list of the held values.

## Details

The check happens at CONSTRUCTION, where the caller can see it, and not
at the fit three layers away. A name the penalty does not carry is an
error naming the ones it does, which is what turns `mcp(x, a = 3)` –
SCAD's shape written on an MCP, whose own is `gamma` – into a message
instead of an argument that lands in `...` and does nothing.

The bounds are OPEN, as they are everywhere in the toolkit: a ridge at
\\\lambda = 0\\ is no penalty at all, and an elastic net at \\\alpha =
0\\ has no kink and is a penalty of another kind.

## See also

[`term_hyper`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md),
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
