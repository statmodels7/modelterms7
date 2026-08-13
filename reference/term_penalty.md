# Penalty of a Term

The penalty attached to the whole of the term's coefficients, or `NULL`
when there is none. The hyperparameters, their bounds and links, and
every derivative in the coefficients and the hyperparameters are the
penalty object's, not the term's.

## Usage

``` r
term_penalty(term, ...)
```

## Arguments

- term:

  An object inheriting from class
  [`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md).

- ...:

  Passed to methods.

## Value

A penalty object, or `NULL`.

## Details

A term whose penalty reaches only part of its parameters returns `NULL`
here and declares that penalty through
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
which names the parameters it covers:
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md)
penalizes the changes and not the linear effect or the break-points, and
the developments of
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`gas`](https://statmodels7.github.io/modelterms7/reference/gas.md)
carry their sub-terms' penalties. Reading a partial penalty here would
say that it covers the block, so the question this generic asks is
answered only where the answer is the whole of it.

## See also

[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
[`term_smooth`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md),
[`edf`](https://statmodels7.github.io/modelterms7/reference/edf.md)

## Examples

``` r
term_penalty(linpar(~x))
#> NULL
```
