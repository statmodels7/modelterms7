# Continuing a Joint Score-Driven Filter Past the Series

Computes the levels that the joint recursion reaches at rows that come
after the observed ones.

## Arguments

- term:

  A built
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- psi:

  The parameters, as for
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- f_past:

  The levels at the observed rows, a matrix with one column per filtered
  parameter.

- s_past:

  The driving scores at the observed rows, a matrix like `f_past`.

- newdata:

  The rows onto which to continue.

- ...:

  Ignored.

- deriv:

  `NULL`, or a list with `df_past` and `ds_past`, lists with one matrix
  per filtered parameter (one row per observed row, `m` columns), and
  `dpsi`, the derivative of the joint parameter vector, one row per
  entry of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  and `m` columns.

## Value

A matrix of continued levels with one row per row of `newdata` and one
column per filtered parameter. With `deriv`, it carries the attribute
`"jacobian"`, a list with one `nrow(newdata)` by `m` matrix per filtered
parameter.

## Details

Past the data each driving score is at its conditional mean, zero, so
the continuation is the deterministic recursion of each level,
\$\$f\_{k,n+h} = \omega_k + \sum_i \Big(\kappa\_{k,i}\\ u\_{k,n+h-i} +
\sum\_{l \in C_k} \kappa\_{k,i,l}\\ u\_{l,n+h-i}\Big) + \sum_j
\xi\_{k,j} f\_{k,n+h-j},\$\$ in which the loadings contribute only while
\\n+h-i\\ is an observed time. A cross loading therefore reads the last
observed scores of the other parameters of its group in the first steps,
after which the levels evolve each on its own. The rules on the new rows
are those of
[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)
for one term: each group continues its own series, the rows of a group
are taken in time order and must come after every observed time of that
group, and the terms must have been built with `time`.

## See also

[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md),
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
