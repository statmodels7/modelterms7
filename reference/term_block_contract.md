# How a Term's Block Moves With Its Coefficients

Returns \\\sum\_{i,j} A\_{ij}\\\partial X\_{ij}/\partial\beta_c\\, one
entry per coefficient of the term: the block's own derivative contracted
against a weight that the caller supplies.

## Usage

``` r
term_block_contract(term, coef = NULL, A, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients, one per coefficient of the term, or `NULL` for the
  coefficients that the term carries.

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
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
does move, and a consumer that differentiates anything built from \\X\\
needs this: a marginal criterion reads \\\log\|K\|\\ with \\K = -\sum_i
w_i \ell\_{ab,i} X_a X_b'\\, so \\\partial K/\partial\beta\\ gains
everything coming from \\\partial X/\partial\beta\\.

Only the contraction is computed. \\\partial X/\partial\beta\\ is \\n
\times m \times m\\ and consumers need only its contractions. For
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) it
is closed form at \\O(nm)\\: writing \\X\_{i,c} = w\_{p}(i) Z_p\[i,c\]\\
with \\w_p = (\partial f/\partial\theta_p)\\h_p'\\, \$\$\frac{\partial
X\_{i,c_1}}{\partial\beta\_{c_2}} =
Z\_{p_1}\[i,c_1\]Z\_{p_2}\[i,c_2\]\Big(f\_{p_1p_2}h\_{p_1}'h\_{p_2}' +
\delta\_{p_1p_2} f\_{p_1} h\_{p_1}''\Big),\$\$ so with \\s\_{p}(i) =
\sum\_{c\in p} A\_{ic}Z\_{p}\[i,c\]\\ the answer is
\\Z\_{p_2}'\sum\_{p_1} q\_{p_1p_2}s\_{p_1}\\, one crossprod per
parameter. The chain rule onto the coefficients, through each
parameter's link and its subformula's design, is applied in this method.

[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
has its own, written from the closed forms. Its derivative is not a
difference that a caller could take: a break-point column is a step
function in its break-point, so the difference quotient diverges as the
step shrinks. The bounded quantities are written out instead. The
truncated line \\(x-\psi)\_+\\ has derivative \\-\mathbf{1}(x\>\psi)\\
in the break-point, and the break-point column
\\-\gamma(x)\mathbf{1}(x\>\psi)\\ has the same indicator as its
derivative in the change and zero almost everywhere in the break-point,
which is the value taken.

A sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
keeps the base method's zeros. Its position is read off a product of the
unknowns, \\\psi = -g/\delta\\, so a column's derivative runs through
that read-off instead of through a development's design, and the weight
\\W = 1/(2\lvert\tilde x-\psi\rvert)\\ that it carries has an unbounded
derivative in the break-point. Its block is a working linearization with
a frozen weight, which is the reason that
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
returns `FALSE` for it. A break-point term with `smoothed =` has a block
that is a true Jacobian, and its contraction is written out from the
smoother's derivatives.

## See also

[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
term_block_contract(b, A = matrix(1, 20, 2))
#> [1] -3.602388  5.798623
```
