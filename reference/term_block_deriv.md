# How a Term's Block Moves Along One Direction

Returns \\\sum_c v_c\\\partial X\_{ij}/\partial\beta_c\\, one entry per
observation and column of the block: the block's own derivative taken
along a direction that the caller supplies.

## Usage

``` r
term_block_deriv(term, coef = NULL, v, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients, one per coefficient of the term, or `NULL` for the
  coefficients that the term carries.

- v:

  A numeric vector as long as the term's coefficients.

- ...:

  Passed to methods.

## Value

A numeric matrix, one row per observation and one column per coefficient
of the term.

## Details

It is the adjoint of
[`term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.md),
and each of the two is computed without the other. That one contracts
over the observations and the columns and returns one entry per
coefficient, as the gradient of a marginal criterion needs; this one
contracts over the coefficients and returns one entry per entry of the
block, as its Hessian needs, because \\\partial K/\partial\beta\\ is
required there in the direction in which the mode moves instead of being
traced.

Both are \\O(nm)\\ and read the same closed form, so this one needs only
the derivatives that the other needs: for
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
writing \\q\_{p_1p_2} = f\_{p_1p_2}h\_{p_1}'h\_{p_2}' +
\delta\_{p_1p_2}f\_{p_1}h\_{p_1}''\\, \$\$\Big(\frac{\partial
X}{\partial\beta}v\Big)\[i, c_1\] = Z\_{p_1}\[i,c_1\]\sum\_{p_2}
q\_{p_1p_2}(i)\\(Z\_{p_2}v\_{p_2})\[i\],\$\$ so the only derivative of
\\f\\ that enters is the second, as in
[`term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.md).

The base method returns zeros, which is exact for a block that does not
move.

## See also

[`term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.md),
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
dim(term_block_deriv(b, v = c(1, 0)))
#> [1] 20  2
```
