# The Second Derivative of the Levinson-Durbin Map

[`gas_levinson`](https://statmodels7.github.io/modelterms7/reference/gas_levinson.md)
with the second derivatives of the coefficients in the partial
autocorrelations propagated as well.

## Usage

``` r
gas_levinson2(pacf)
```

## Arguments

- pacf:

  A numeric vector of partial autocorrelations in \\(-1, 1)\\.

## Value

A list with `phi`, `jacobian` and `hessian`, the last a list of one
symmetric matrix per coefficient.

## Details

The recursion is \$\$\phi^{(k)}\_k = \rho_k, \qquad \phi^{(k)}\_i =
\phi^{(k-1)}\_i - \rho_k\phi^{(k-1)}\_{k-i},\$\$ which is bilinear:
\\\rho_k\\ multiplies quantities that do not depend on it.
Differentiating twice therefore adds no new kind of term, only the two
places the product rule puts the first derivative, \$\$H^{(k)}\_i =
H^{(k-1)}\_i - \rho_k H^{(k-1)}\_{k-i} - e_k
\left(J^{(k-1)}\_{k-i}\right)^{\\\top} - J^{(k-1)}\_{k-i}
e_k^{\top},\$\$ and the last coefficient's second derivative is zero at
every order, it being \\\rho_k\\ itself.

It is wanted because the observed information of a model carrying this
term needs the second derivative of the predictor in the term's own
parameters, and the persistence reaches the predictor through this map.

## See also

[`gas_levinson`](https://statmodels7.github.io/modelterms7/reference/gas_levinson.md)
