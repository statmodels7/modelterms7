# The Fourth Derivative of a Product, in Two Directions

Computes the fourth derivative of a product as a sum of sixteen terms,
one for each way of distributing four differentiations between two
factors, where two of the four slots are contracted against the
directions.

## Usage

``` r
.gas_prod4(A, B)
```

## Arguments

- A, B:

  The two factors' jets, as
  [`.gas_jet()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_jet.md)
  builds them.

## Value

A matrix, the contracted fourth derivative of the product.

## Details

The recursion forms products of three kinds: a loading times a lagged
score, \\\kappa_i s\_{t-i}\\, a persistence times a lagged level,
\\\xi_j f\_{t-j}\\, and the starting level, which is the fixed point
\\f_0 = \omega + Sf_0\\. The rule is therefore written once here. Let
\\S\\ be the set of slots that go to the first factor. Once \\m\\ is
contracted against \\v\\ and \\n\\ against \\w\\, the sixteen subsets of
\\\\m,n,b,c\\\\ give the terms that this function returns. The sum is
symmetric in the two factors, because the swap maps each term onto
another term.

## See also

[`.gas_jet()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_jet.md)
