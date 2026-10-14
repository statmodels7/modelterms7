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
two directions in which the mode moves, and the persistence reaches the
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

Both third-order contractions \\T\[v\]\\ and \\T\[w\]\\ are carried: the
terms in \\v_k\\ and \\w_k\\ use one each, and the last two terms use
\\T\[v\]\\ against \\w\\, so a recursion holding a single contraction
cannot reach this order.

Each coefficient \\\xi^{(k)}\_i\\ is a polynomial of degree at most
\\k\\ in the first \\k\\ partial autocorrelations, multilinear in them,
so the result is identically zero for \\q \le 3\\. The first coefficient
carrying a monomial of degree four is \\\xi^{(4)}\_2\\, which needs \\q
= 4\\. A comparison with a reference therefore needs \\q \ge 4\\.

## See also

[`gas_levinson3()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson3.md),
[`gas_levinson2()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson2.md)
