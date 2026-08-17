# How a Term's Block Moves With Its Coefficients

\\\sum\_{i,j} A\_{ij}\\\partial X\_{ij}/\partial\beta_c\\, one entry per
coefficient of the term: the block's own derivative contracted against a
weight the caller supplies.

## Usage

``` r
term_block_contract(term, coef = NULL, A, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients, or `NULL` for the ones it carries.

- A:

  A numeric matrix, one row per observation and one column per
  coefficient of the term.

- ...:

  Passed to methods.

## Value

A numeric vector as long as the term's coefficients.

## Details

A term whose block is a fixed design does not move, and the base method
returns zeros. A term registering
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
does move, and a consumer that differentiates anything built from \\X\\
needs this: a marginal criterion reads \\\log\|K\|\\ with \\K = -\sum_i
w_i \ell\_{ab,i} X_a X_b'\\, so \\\partial K/\partial\beta\\ gains
everything coming from \\\partial X/\partial\beta\\.

**The contraction and not the derivative**, because \\\partial
X/\partial\beta\\ is \\n \times m \times m\\ and nothing needs it whole.
For [`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md) it
is closed form at \\O(nm)\\: writing \\X\_{i,c} = w\_{p}(i) Z_p\[i,c\]\\
with \\w_p = (\partial f/\partial\theta_p)\\h_p'\\, \$\$\frac{\partial
X\_{i,c_1}}{\partial\beta\_{c_2}} =
Z\_{p_1}\[i,c_1\]Z\_{p_2}\[i,c_2\]\Big(f\_{p_1p_2}h\_{p_1}'h\_{p_2}' +
\delta\_{p_1p_2} f\_{p_1} h\_{p_1}''\Big),\$\$ so with \\s\_{p}(i) =
\sum\_{c\in p} A\_{ic}Z\_{p}\[i,c\]\\ the answer is
\\Z\_{p_2}'\sum\_{p_1} q\_{p_1p_2}s\_{p_1}\\, one crossprod per
parameter. The chain rule onto the coefficients – each parameter's link
and its subformula's design – stays here, which is the only place that
knows them.

⚠️ The base method is right for a FIXED block and is what
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
inherit, whose blocks do move. Theirs is not a difference a caller could
take either: a break-point column is a step function in its break-point,
so the quotient diverges as the step shrinks – measured at h, h/4 and
h/16, 3.6e4, 1.4e5 and 5.8e5, against
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md)'s
0.6038 throughout. They need their own method, written from the closed
forms, which is not this release.

## See also

[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`nl_fderiv`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
term_block_contract(b, A = matrix(1, 20, 2))
#> [1] -3.602388  5.798623
```
