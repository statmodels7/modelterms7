# How a Term's Block Moves Along Three Directions

Returns \\\sum\_{c,d,e} v_c u_d s_e\\\partial^3 X\_{ij}/\partial\beta_c
\partial\beta_d\partial\beta_e\\, one entry per observation and column
of the block: the third derivative of a term's block along three
directions. The Jeffreys penalty of a block inside the term's subformula
reads it in the derivative of the outer criterion.

## Usage

``` r
term_block_deriv3(term, coef = NULL, v, u, s, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients, one per coefficient of the term, or `NULL` for the
  coefficients that the term carries.

- v, u, s:

  Numeric vectors as long as the term's coefficients.

- ...:

  Passed to methods.

## Value

A numeric matrix, one row per observation and one column per coefficient
of the term.

## Details

The base method returns zeros, which is exact for a block that does not
move. For
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
with \\\theta_p = h_p(Z_p\beta_p)\\ and the block's column \\c\\ in
parameter \\p_1\\, the entry is the fourth derivative of the predictor
along \\(e_c, v, u, s)\\, the sum over the fifteen partitions of the
four directions, each block of a partition carrying the derivative of
one link and each partition the derivative of \\f\\ of its number of
blocks. Writing \\a_q, b_q, c_q\\ for the three directions carried onto
parameter \\q\\ (\\Z_q v_q\\ and so on), \\h^{(k)}\_q\\ for the
derivatives of its inverse link and \\f\_{q\dots}\\ for those of \\f\\,
\$\$\begin{aligned} X^{(3)}\[t,c\] =
Z\_{p_1}\[t,c\]\Big\[&f\_{p_1}h^{(4)}\_{p_1}a\_{p_1}b\_{p_1}c\_{p_1} +
h^{(3)}\_{p_1}\big(a\_{p_1}b\_{p_1}S(c) + a\_{p_1}c\_{p_1}S(b) +
b\_{p_1}c\_{p_1}S(a)\big)\\ &+ h''\_{p_1}\big(a\_{p_1}(Q(b,c) +
R(b,c)) + b\_{p_1}(Q(a,c) + R(a,c)) + c\_{p_1}(Q(a,b) + R(a,b))\big)\\
&+ h'\_{p_1}\big(\textstyle\sum_q f\_{p_1q}h^{(3)}\_qa_qb_qc_q +
P(ab,c) + P(ac,b) + P(bc,a) +
\sum\_{q,r,s}f\_{p_1qrs}h'\_qa_q\\h'\_rb_r\\h'\_sc_s\big)\Big\],
\end{aligned}\$\$ with \\S(x) = \sum_q f\_{p_1q}h'\_qx_q\\, \\Q(x,y) =
\sum_q f\_{p_1q}h''\_qx_qy_q\\, \\R(x,y) =
\sum\_{q,r}f\_{p_1qr}h'\_qx_q\\h'\_ry_r\\ and \\P(xy,z) =
\sum\_{q,r}f\_{p_1qr}h''\_qx_qy_q\\h'\_rz_r\\. A parameter whose
development moves with the coefficients has a design of its own that
moves, and its third derivative is not written: the method signals an
error there.

## See also

[`term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.md),
[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md)

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
dim(term_block_deriv3(b, v = c(1, 0), u = c(0, 1), s = c(1, 1)))
#> [1] 20  2
```
