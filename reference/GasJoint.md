# A Joint Score-Driven Filter

The S7 class of several
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
terms that sit in the equations of different distribution parameters and
are filtered as one recursion. Build one with
[`gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.md).

## Usage

``` r
GasJoint(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  ids = character(0),
  blueprint = list(),
  terms = list(),
  params = character(0),
  offsets = integer(0),
  partners = list(),
  keys = character(0)
)
```

## Arguments

- label:

  A character string prefixed to the term's coefficient names when
  non-empty, and used as the title of
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and the tag
  of [`print()`](https://rdrr.io/r/base/print.html). The empty string
  `""` means no label.

- hyper:

  The hyperparameters of the term's penalty that the caller held, as a
  named list keyed by the penalty's names. Empty, the default, means
  every one of them is estimated. See
  [`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md).

- grid:

  How many values a path visits for each of the term's hyperparameters,
  as a named list of single whole numbers. Empty, the default, leaves
  the number to the fitting layer. See
  [`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md).

- values:

  The values that a path visits, for each hyperparameter that the caller
  wrote out, as a named list of numeric vectors. Empty, the default, has
  the path build them. See
  [`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md).

- min_ratio:

  How far down the path over the size of the kink reaches, as a fraction
  of the value that empties the block: one number in \\(0, 1)\\, or
  `numeric(0)` for the fitting layer's own. See
  [`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md).

- search:

  How the term's own hyperparameters are covered when it has more than
  one carrying a kink: `"grid"` for every combination of them,
  `"cyclic"` for one at a time, or `character(0)` for the default. See
  [`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md).

- ids:

  Which of the term's hyperparameters are shared with those of other
  terms, and under what label: a character vector named by the term's
  own hyperparameters, or `character(0)` for none. The terms carrying
  the same label for the same hyperparameter estimate one value. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md).

- blueprint:

  The shared ordering: `order` (the rows of each group, in time order)
  and `n`. The property is inherited from
  [`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md).

- terms:

  A named list of built
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md)
  objects, one per filtered parameter, named by the parameter.

- params:

  The names of the filtered parameters, in the order of `terms`.

- offsets:

  The position before the first parameter of each term in the joint
  parameter vector.

- partners:

  A list with one integer vector per term: the positions in `params` of
  the parameters whose scores the term's level loads besides its own.

- keys:

  The names under which a fitting layer knew each term before joining
  them, named by the parameter, or an empty vector.

## Value

An S7 object of class `GasJoint`, inheriting from
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md),
so that a fitting layer treats the joint filter as one structural term
whose parameters are those of its terms, each prefixed by the
distribution parameter of its term.

## Details

A score-driven level in the equation of one parameter is driven by the
score of that parameter, which is evaluated at the current values of all
the others. Where two parameters both carry a level, each score depends
on the other level, so the two recursions are one: with \\G\\ the
filtered parameters, \$\$f\_{k,t} = \omega_k + \sum\_{i=1}^{p_k}
\Big(\kappa\_{k,i}\\ u\_{k,t-i} + \sum\_{l \in C_k} \kappa\_{k,i,l}\\
u\_{l,t-i}\Big) + \sum\_{j=1}^{q_k} \xi\_{k,j} f\_{k,t-j}, \qquad k \in
G,\$\$ where \\u_t\\ is the vector of driving scores of the parameters
in \\G\\, read at the predictor of every equation at time \\t\\, and
\\C_k\\ is the set of the other parameters of the group of `cross` that
contains \\k\\ (empty when the term carries no label). The persistence
of each level is its own, so the autoregressive part is stationary when
each equation's is.

The driving vector is the caller's: a fitting layer supplies \\u_t =
\mathcal{I}\_{GG,t}^{-d}\\\partial\ell_t/\partial\eta\_{G,t}\\ with the
expected information of the filtered parameters and its derivative in
the predictors.

## See also

[`gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.md),
the constructor;
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
for one term.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = rnorm(40))
tm <- list(mu = term_build(gas(p = 1, q = 1, time = t), dd),
           sigma = term_build(gas(p = 1, q = 1, time = t), dd))
jt <- gas_joint(tm)
S7::S7_inherits(jt, GasJoint)
#> [1] TRUE
jt@params
#> [1] "mu"    "sigma"
```
