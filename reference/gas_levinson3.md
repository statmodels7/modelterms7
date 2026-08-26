# The Third Derivative of the Levinson-Durbin Map, in One Direction

[`gas_levinson2()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson2.md)'s
second derivatives differentiated once more and contracted against a
single direction, one matrix per coefficient.

## Usage

``` r
gas_levinson3(pacf, w)
```

## Arguments

- pacf:

  A numeric vector of partial autocorrelations in \\(-1, 1)\\.

- w:

  The direction to contract against, as long as `pacf`.

## Value

A list of one `q` by `q` matrix per coefficient.

## Details

The exact gradient of a marginal criterion over a penalty on this term's
own parameters needs the third derivative of the predictor, and the
persistence reaches the predictor through this map. It is needed only
contracted: the criterion asks for \\\mathrm{tr}(M\\\partial K/\partial
u\[v\])\\, a derivative along the single direction the penalized mode
moves in, so what is propagated is a matrix per coefficient and never a
three-index array.

Differentiating the hessian recursion of
[`gas_levinson2()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson2.md)
once more along \\w\\ adds no new kind of term, the map being bilinear:
\$\$T^{(k)}\_i = T^{(k-1)}\_i - \rho_k T^{(k-1)}\_{k-i} - w_k
H^{(k-1)}\_{k-i} - e_k\left(H^{(k-1)}\_{k-i}w\right)^{\\\top} -
\left(H^{(k-1)}\_{k-i}w\right)e_k^{\top},\$\$ and the last coefficient's
third derivative is zero at every order, it being \\\rho_k\\ itself.

The map is multilinear of degree \\k\\ in the first \\k\\ partial
autocorrelations, so the result is identically zero for \\q \le 2\\: at
\\q = 2\\ the only non-trivial coefficient is \\\phi_1 =
\rho_1(1-\rho_2)\\, which is bilinear. A check of this function that
stops at \\q = 2\\ compares zero with zero and asserts nothing.

## See also

[`gas_levinson2()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson2.md)
