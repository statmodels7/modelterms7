# Every Penalty That a Term Carries

Returns the penalties that a term declares: a list of entries, where
each entry names a subset of the term's own parameters, the penalty over
them, and the hyperparameters that the caller held for that penalty. A
fitting layer enumerates penalties through this list, and
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
and
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
read it.

## Usage

``` r
term_penalties(term, ...)
```

## Arguments

- term:

  A built term. An unbuilt one is accepted, and the result lists the
  entries that it already declares, usually none.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

An unnamed list, possibly empty, one element per declared penalty. Every
element is a list with the same nine fields; an empty value stands for a
field that is not set:

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

  the hyperparameters that the caller held, a named list, empty when all
  of them are estimated.

- `n_values`, `values`, `min_ratio`, `search`:

  what the caller specified about the path over this penalty's
  hyperparameters, taken from the term's properties `grid`, `values`,
  `min_ratio` and `search`.

- `ids`:

  the sharing labels of this penalty's hyperparameters, a named
  character vector, empty where the penalty shares none of its
  hyperparameters.

An entry that a term carries up from one of its sub-terms keeps the
values and the search that the sub-term holds.

## Two ways in which it generalizes [`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)

A term may carry more than one penalty, over different parameters of its
own. A panel model with a population value and a departure per group
leaves the population value free and shrinks the departures, which is
one penalty over part of the parameters and no penalty over the rest. An
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) term
that develops two of its parameters by two different penalized sub-terms
carries two entries.

The parameters need not be coefficients of a design block. The
persistence of a score-driven term, the nonlinear parameters of
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
the break-point of
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
are parameters of the term and of nothing else, and a penalty needs from
them only a vector of numbers and their positions. For a structural term
`index` gives positions in
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md);
for an additive one, columns of the block.

## The base method, and the name of an entry

The method on
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
derives its result from
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
so a term carrying one penalty over its whole block needs no method of
its own. Its single entry is named `""`, meaning the whole term.

A name is unique within the term and is not the term's own name. Two
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
terms in one formula are two terms with their own hyperparameters, told
apart by their names in the formula;
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
records a hyperparameter under the key `term` for a term with one entry
and `term::entry` for a term with several.

The list itself is unnamed: the key of an entry is its `name` field,
read as `e$name`, and `names(entries)` is `NULL`.

## A specification carries no entries

The penalty is attached at
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
so `term_penalties()` on an unbuilt penalized term is an empty list.
Consequently
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md),
which reads this list, returns `TRUE` for `lasso(~ x)` and `FALSE` once
the term is built.

## See also

[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
for the single-penalty case,
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
for how a term's columns divide among its own parameters,
[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)
for only the held values, and
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what the entries cost.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(50, 0, 10)), g = factor(rep(c("a", "b"), 25)))

# A penalty over the whole block gives one entry, named with the empty
# string.
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
