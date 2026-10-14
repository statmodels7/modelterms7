# The Fourth Derivative of a Product, in Two Directions

Sixteen terms, one per way of distributing four differentiations between
two factors, with two of the four slots contracted against the
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
\\f_0 = \omega + Sf_0\\. The rule is therefore written once here.
Writing \\S\\ for the set of slots that go to the first factor, the
sixteen subsets of \\\\m,n,b,c\\\\ give, after contracting \\m\\ against
\\v\\ and \\n\\ against \\w\\, the terms this returns. The sum is
symmetric in the two factors, each term mapping onto another under the
swap.

## See also

[`.gas_jet()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_jet.md)
