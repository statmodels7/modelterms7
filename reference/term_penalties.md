# Every Penalty a Term Carries

What a term declares it wants penalized: a list of entries, each naming
a subset of the term's own parameters and the penalty over them.

## Usage

``` r
term_penalties(term, ...)
```

## Arguments

- term:

  A built term.

- ...:

  Passed to methods.

## Value

A list, possibly empty. Each entry has `name` (a label unique WITHIN the
term, empty for a penalty over the whole of it), `index` (positions
among the term's parameters) and `penalty` (a penalties7 object). The
name is not the term's: two
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
terms in one formula are two terms with their own hyperparameters, and
it is the caller that knows what it called each one.

## Details

[`term_penalty`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
answers for the common case, one penalty over the whole of a term's
design block, and this generalizes it in two directions a model layer
needs.

A term may carry **more than one** penalty, over different parameters of
its own. A panel model with a population value and a departure per group
wants the population value free and the departures shrunk, which is one
penalty over part of the parameters and none over the rest.

The parameters need **not be coefficients of a design block**. The
persistence of a score-driven term, the nonlinear parameters of
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md), the
break-point of
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md) are
parameters of the term and nothing else, and everything a penalty needs
from them is a vector of numbers and their positions.

The base method answers from
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
so a term that carries one penalty over its whole block – every term
shipped here – needs no method of its own and behaves exactly as before.
Its single entry is named with the empty string, meaning the whole term,
so a caller that keys the hyperparameters by term name keys them exactly
as it did.

## See also

[`term_penalty`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
[`term_npar`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)

## Examples

``` r
term_penalties(term_build(ridge(~x), data.frame(x = rnorm(20))))
#> [[1]]
#> [[1]]$name
#> [1] ""
#> 
#> [[1]]$index
#> [1] 1
#> 
#> [[1]]$penalty
#> quadratic penalty on 1 coefficient(s) through 1 row(s); theta: lambda
#> 
#> [[1]]$fixed
#> list()
#> 
#> 
term_penalties(term_build(linpar(~x), data.frame(x = rnorm(20))))
#> list()
```
