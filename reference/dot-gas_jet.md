# A Quantity's Derivatives and Their Contractions

Collects a scalar quantity's derivatives through the fourth order,
together with the contractions of the lower ones against the two
directions, in the shape
[`.gas_prod4()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_prod4.md)
reads.

## Usage

``` r
.gas_jet(v0, d1, d2, d3v, d3w, d4, vks)
```

## Arguments

- v0:

  The value.

- d1, d2:

  The first and second derivatives.

- d3v, d3w:

  The third derivative contracted against each direction.

- d4:

  The fourth derivative contracted against both.

- vks:

  The two directions, on the active set.

## Value

A list with `v`, `d1`, `d2`, `d3v`, `d3w`, `d4` and the contractions
`a`, `b`, `d2v`, `d2w`, `d2vw` and `d3vw`.

## Details

The fourth derivative of a product needs, of each factor, its value, its
first three derivatives and its fourth contracted against both
directions, and also the contractions the product rule pairs against the
other factor's. Computing them once per factor rather than at each of
the sixteen terms is what keeps the recursion's cost the same shape as
the third order's.

## See also

[`.gas_prod4()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_prod4.md)
