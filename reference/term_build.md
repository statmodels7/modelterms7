# Build a Term on Data

Turns a term specification into a built term: the design block is
computed from the data, the coefficient names are assigned, and the
blueprint that reproduces the mapping on new data is recorded. The
returned object is a copy of the specification with those properties
filled; the specification itself is unchanged.

## Usage

``` r
term_build(term, data, ...)
```

## Arguments

- term:

  An object inheriting from class
  [`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md).

- data:

  A data frame.

- ...:

  Passed to methods.

## Value

A built term of the same class as `term`.

## Details

An additive term contributes to the linear predictor through a design
block and, when it is penalized, a penalty on the coefficients of that
block:

\$\$\eta = \sum\_{t} X_t \beta_t, \qquad \text{penalized objective}
\quad -\ell(\beta) + \sum\_{t} \rho_t(\beta_t; \theta_t),\$\$

and building the term is what produces \\X_t\\ from the data and
attaches \\\rho_t\\.
[`term_matrix`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
reads the block,
[`term_penalty`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
the penalty and
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reproduces \\X_t\\ on new rows through the blueprint. A structural term
is the exception: its contribution cannot be written as a block of
columns, and it reports itself through
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
or
[`term_loglik`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
instead.

## Examples

``` r
built <- term_build(linpar(~x), data.frame(x = 1:4))
term_matrix(built)
#>   (Intercept) x
#> 1           1 1
#> 2           1 2
#> 3           1 3
#> 4           1 4
```
