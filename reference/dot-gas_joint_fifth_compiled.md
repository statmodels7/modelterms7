# The Fifth Order of a Joint Filter, Compiled

Builds the groups that `gas_fifth_joint_cpp()` reads and runs it, for a
joint filter whose members carry subformulas or not.

## Usage

``` r
.gas_joint_fifth_compiled(term, psi, g, seed, blocks_data, dirs, eta = NULL)
```

## Arguments

- term:

  A built term.

- psi:

  The term's parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The weights, one column per filtered parameter.

- seed:

  A list of one matrix of static rows per equation.

- blocks_data:

  The model's derivative arrays, as set out above.

- dirs:

  The three directions.

- eta:

  The static part of the predictor.

## Value

As
[`term_fifth.GasJoint()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.GasJoint.md).
