# Which of a Term's Hyperparameters Are Shared

The sharing labels a term carries. Each element of the answer is a named
character vector whose names are the hyperparameters of one penalty
entry and whose values are the labels. The terms carrying the same label
for the same hyperparameter estimate **one** value.

## Usage

``` r
term_ids(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A list, one named character vector per penalty entry that shares
anything, keyed by that entry's name. Empty where the term shares
nothing.

## What sharing is, and what it is not

Sharing does not merge the penalties. They stay two objects with the
same number in them, so each keeps its own block, its own coefficients
and its own effective degrees of freedom; what is shared is the value
that a criterion estimates. Two smooths with one label are one smoothing
parameter over two curves, which is mgcv's `id`, and two lassos with one
label are one threshold swept along one path.

## Comparability of the penalties

The comparability of two penalties is not checked, deliberately. A
shared value multiplies each penalty as it is, so it has the usual
meaning for penalties of the same shape over covariates of comparable
scale, and the caller is responsible for the comparison in other cases.
The intended case is two smooths of the same basis and dimension.

## Where it may be written

On any term that carries a penalty, including one written inside a
subformula: a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
developing its slope change and its break-point over two ridges may give
the two one label and estimate one shrinkage for them, and a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
developing its slope change and its level change may do the same.

The shape is that of
[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md):
one element of the answer for each entry of
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
that shares anything, keyed by the name of the entry, so the two
functions read alike. A built term answers from its entries, which
allows a structural term to report what the terms of its subformulas
share without a method of its own, and an unbuilt term answers from its
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
