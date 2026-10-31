# The Fifth Derivative of the Levinson-Durbin Map, in Three Directions

Differentiates the fourth derivatives of
[`gas_levinson4()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson4.md)
once more and contracts them against a third direction, which gives one
matrix per coefficient.

## Usage

``` r
gas_levinson5(pacf, v, w, x)
```

## Arguments

- pacf:

  A numeric vector of partial autocorrelations in \\(-1, 1)\\.

- v, w, x:

  The three directions, each as long as `pacf`.

## Value

A list of one `q` by `q` matrix per coefficient.

## Details

The map is bilinear at every step of the recursion, so differentiating
it once more adds no new kind of term. Contracting three of the five
slots of the full fifth derivative against \\v\\, \\w\\ and \\x\\ leaves
\$\$F^{(k)}\_i = F^{(k-1)}\_i - \rho_k F^{(k-1)}\_{k-i} - v_k
Q^{(k-1)}\_{k-i}\[w,x\] - w_k Q^{(k-1)}\_{k-i}\[v,x\] - x_k
Q^{(k-1)}\_{k-i}\[v,w\] - e_k\big(Q^{(k-1)}\_{k-i}\[v,w\]x\big)^{\top} -
\big(Q^{(k-1)}\_{k-i}\[v,w\]x\big)e_k^{\top},\$\$ where
\\Q\[\cdot,\cdot\]\\ are the fourth-order contractions of
[`gas_levinson4()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson4.md),
all three of which are carried, with the third-order contractions along
each direction that they need. The result is symmetric in the three
directions, and identically zero for \\q \le 4\\.

## See also

[`gas_levinson4()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson4.md)
