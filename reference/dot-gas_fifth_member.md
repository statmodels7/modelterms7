# The Chart of One Member for the Compiled Fifth Order

Returns, for a built score-driven term, the values that the compiled
kernels read per row: the level, the loadings, the partial
autocorrelations and the cross loadings, each as design rows over the
unknowns with its link's derivatives of orders zero to five.
`.gas_fifth_rows()` restricts them to a group's rows and active
unknowns.

## Usage

``` r
.gas_fifth_rows(mem, rows, act)

.gas_fifth_member(tm, psiv, cols, dirs)
```

## Arguments

- mem:

  The result of `.gas_fifth_member()`.

- rows, act:

  A group's rows and active unknowns.

- tm:

  A built score-driven term.

- psiv:

  Its parameters in the order of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- cols:

  The columns of those parameters among the unknowns.

- dirs:

  The three directions over the unknowns.

## Value

A list with `vals` (named by the layout's slots), `p`, `q`, `cross` (the
cross loadings' names), `varying_b` (whether a partial autocorrelation
carries a subformula) and `active`, the term's coordinates active in
each group.
