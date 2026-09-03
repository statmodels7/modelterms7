# Which of a Term's Hyperparameters Are Shared

The sharing labels a term carries: a named character vector whose names
are its own hyperparameters and whose values are the labels. The terms
carrying the same label for the same hyperparameter estimate **one**
value.

## Usage

``` r
term_ids(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A list, one named character vector per penalty entry that shares
anything, keyed by that entry's name. Empty where the term shares
nothing.

## What sharing is, and what it is not

It does not merge the penalties. They stay two objects with the same
number in them, so each keeps its own block, its own coefficients and
its own effective degrees of freedom; what is shared is the value a
criterion estimates. Two smooths with one label are one smoothing
parameter over two curves, which is mgcv's `id`, and two lassos with one
label are one threshold swept along one path.

## Care

Nothing checks whether two penalties are comparable, deliberately. A
shared value multiplies each penalty as it is, so it means what one
expects between penalties of the same shape over covariates of
comparable scale, and something the caller had better intend otherwise.
Two smooths of the same basis and dimension are the case it is for.

## Where it may be written

On any term that carries a penalty, including one written inside a
subformula: a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
developing its slope change, its level change and its break-point over
three ridges may give the three one label and estimate one shrinkage for
them.

The shape is
[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)'s,
one entry of the answer per entry of
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
so the two reporters read alike: a built term answers from its entries,
which is what lets a structural term report what its subformulas' terms
share without a method of its own, and an unbuilt one answers from its
own property.

## See also

[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)
for holding a hyperparameter instead of sharing it,
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the entries the labels belong to.

## Examples

``` r
# One hyperparameter, so the label needs no name.
term_ids(s(x, id = "L"))
#> [[1]]
#> lambda 
#>    "L" 
#> 

# Several, so it does: this shares the mixing and leaves the strength free.
term_ids(enet(~ x, id = c(alpha = "A")))
#> [[1]]
#> alpha 
#>   "A" 
#> 

# A term sharing nothing answers with an empty list.
term_ids(s(x))
#> list()
term_ids(linpar(~ x))
#> list()
```
