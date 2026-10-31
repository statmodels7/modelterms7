# The Autoregressive Coefficients' Jets from Rank-One Partial Autocorrelations

Returns the jets of the coefficients \\\xi_j(\rho)\\ of the
Levinson-Durbin map, where each partial autocorrelation is \\\rho_k =
h_k(z_k^\top u)\\.

## Usage

``` r
.gas_map_jets(K, zs, dirs)
```

## Arguments

- K:

  A `q` by 6 matrix: the inverse links and their derivatives of orders
  one to five, per partial autocorrelation.

- zs:

  A `q` by `m` matrix of the rows \\z_k\\.

- dirs:

  The directions.

## Value

A list of `q` jets.

## Details

The derivative of order \\n\\ of \\\rho_k\\ is \\h_k^{(n)} z_k^{\otimes
n}\\, so Faa di Bruno's formula over the set partitions of the slots
reads, for each block, a vector over the partial autocorrelations (a
block without a free slot), a \\q \times m\\ matrix (a block with one
free slot) or a weighted sum of \\z_k z_k^\top\\ (a block with both).
The map's derivatives come from
[`.gas_levinson_contract()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_levinson_contract.md).
With \\z_k\\ a unit vector this is
[`.gas_chart_contract()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_contract.md)
lifted to the unknowns.
