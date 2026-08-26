# Every Penalty a Term Carries

Returns what a term declares it wants penalized: a list of entries, each
naming a subset of the term's own parameters, the penalty over them, and
whatever the caller fixed about that penalty's hyperparameters. This is
the enumeration a fitting layer runs over, and it is what
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
and
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
read.

## Usage

``` r
term_penalties(term, ...)
```

## Arguments

- term:

  A built term. An unbuilt one is accepted and reports what it has,
  which is usually nothing.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

An unnamed list, possibly empty, one element per declared penalty. Each
element is a list with

- `name`:

  a character label unique within the term, `""` for a penalty covering
  the whole of it.

- `index`:

  integer positions among the term's own parameters: columns of the
  block for an additive term, positions in
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  for a structural one.

- `penalty`:

  a penalties7 penalty over exactly those parameters, so
  `penalty@n_coef` equals `length(index)`.

- `fixed`:

  the hyperparameters the caller held, a named list, empty when all of
  them are estimated.

- `n_values`, `values`, `min_ratio`, `search`:

  what the caller said about the path over this penalty's
  hyperparameters, from the term's properties of the same names.

## Two ways it generalizes [`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)

**A term may carry more than one penalty**, over different parameters of
its own. A panel model with a population value and a departure per group
wants the population value free and the departures shrunk, which is one
penalty over part of the parameters and none over the rest.
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
developing two of its parameters by two different penalized sub-terms
carries two entries.

**The parameters need not be coefficients of a design block.** The
persistence of a score-driven term, the nonlinear parameters of
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
the break-point of
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
are parameters of the term and of nothing else, and all a penalty needs
from them is a vector of numbers and their positions. For a structural
term `index` gives positions in
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md);
for an additive one, columns of the block.

## The base method, and the name of an entry

The method on
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
answers from
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
so a term carrying one penalty over its whole block needs no method of
its own. Its single entry is named `""`, meaning the whole term.

A name is unique **within** the term and is not the term's own name. Two
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
terms in one formula are two terms with their own hyperparameters, and
it is the caller who knows what it called each of them;
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
composes a key as `term` or `term::entry`.

The list itself is **unnamed**: read `e$name`, not `names(entries)`.

## A specification carries no entries

The penalty is attached at
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
so `term_penalties()` on an unbuilt penalized term is an empty list.
That is why
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md),
which reads this, answers `TRUE` for `lasso(~ x)` and `FALSE` once it is
built.

## See also

[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
for the single-penalty case,
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
for how a term's columns divide among its own parameters,
[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)
for the held values alone, and
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what the entries cost.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(50, 0, 10)), g = factor(rep(c("a", "b"), 25)))

# One penalty over the whole block: one entry, named with the empty string.
e <- term_penalties(term_build(ridge(~ x), d))
length(e)
#> [1] 1
str(e[[1]][c("name", "index")])
#> List of 2
#>  $ name : chr ""
#>  $ index: int 1
e[[1]]$penalty
#> quadratic penalty on 1 coefficient(s) through 1 row(s); theta: lambda

# An unpenalized term declares nothing.
term_penalties(term_build(linpar(~ x), d))
#> list()

# A penalty over part of a term's parameters: nl() with one of its two
# parameters developed by a lasso, so the entry covers columns 1 and 2
# of a block of three.
nb <- term_build(nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g),
                    start = list(r = 1.3)), d)
ent <- term_penalties(nb)
vapply(ent, function(z) z$name, character(1))
#> [1] "a::lasso(~g)"
ent[[1]]$index
#> [1] 1 2
term_npar(nb)
#> [1] 3

# The list is unnamed: the key is the entry's own field.
names(ent)
#> NULL

# Unbuilt, there is nothing to report yet.
term_penalties(lasso(~ x))
#> list()
```
