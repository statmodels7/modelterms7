# The Fourth Derivative of the Levinson-Durbin Map, in Two Directions

[`gas_levinson3()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson3.md)'s
third derivatives differentiated once more and contracted against a
second direction, one matrix per coefficient.

## Usage

``` r
gas_levinson4(pacf, v, w)
```

## Arguments

- pacf:

  A numeric vector of partial autocorrelations in \\(-1, 1)\\.

- v, w:

  The two directions to contract against, each as long as `pacf`. The
  result is symmetric in the two.

## Value

A list of one `q` by `q` matrix per coefficient.

## Details

The second derivative of a marginal criterion in a pair of
hyperparameters reads the fourth derivative of the predictor along the
two directions the mode moves in, and the persistence reaches the
predictor through this map. As at the third order it is needed only
contracted, so what is propagated is a matrix per coefficient and never
a four-index array.

Differentiating the third-order recursion of
[`gas_levinson3()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson3.md)
once more adds no new kind of term, the map being bilinear at every
step. The full fourth derivative obeys \$\$Q^{(k)}\_i = Q^{(k-1)}\_i -
\rho_k Q^{(k-1)}\_{k-i} - \sum\_{\text{4 slots}} \delta\_{\cdot
k}T^{(k-1)}\_{k-i}\[\text{the rest}\],\$\$ and contracting two of those
slots against \\v\\ and \\w\\ leaves \$\$Q^{(k)}\_i\[v,w\] =
Q^{(k-1)}\_i\[v,w\] - \rho_k Q^{(k-1)}\_{k-i}\[v,w\] - v_k
T^{(k-1)}\_{k-i}\[w\] - w_k T^{(k-1)}\_{k-i}\[v\] -
e_k\big(T^{(k-1)}\_{k-i}\[v\]w\big)^{\\\top} -
\big(T^{(k-1)}\_{k-i}\[v\]w\big)e_k^{\top}.\$\$

Both third derivatives are carried, not one: the last two terms couple
\\T\[v\]\\ with \\w\\ and \\T\[w\]\\ with \\v\\, so a recursion holding
a single contraction cannot reach this order.

The map is multilinear of degree \\k\\ in the first \\k\\ partial
autocorrelations, so the result is identically zero for \\q \le 3\\: the
first coefficient carrying a monomial of degree four is
\\\phi^{(4)}\_1\\, which needs \\q = 4\\. A check of this function that
stops at \\q = 3\\ compares zero with zero and asserts nothing.

## See also

[`gas_levinson3()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson3.md),
[`gas_levinson2()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson2.md)
