# The Second Derivative of a Design Block in the Coefficients

Returns \\\sum\_{q,r} v_q u_r\\\partial^2
X\_{ij}/\partial\beta_q\partial\beta_r\\, one entry per observation and
column of the block: the block's second derivative contracted in two
directions that the caller supplies.

## Usage

``` r
term_block_deriv2(term, coef = NULL, v, u, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients, one per coefficient of the term, or `NULL` for the
  coefficients that the term carries.

- v, u:

  Numeric vectors as long as the term's coefficients.

- ...:

  Passed to methods.

## Value

A numeric matrix, one row per observation and one column per coefficient
of the term.

## Details

It stands one order above
[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md),
and it is contracted, never returned as an array, for the same reason as
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
in the structural branch: the full object has \\m\\ coefficient indices
twice over beside the \\n\\ rows and \\m\\ columns of the block, and
only its contraction in the two directions in which the penalized mode
moves is ever read.

The quantity enters the Hessian of a marginal criterion and nothing
else. It is absent from the criterion's value, from its gradient and
from the fit, so the coefficients, the log-likelihood, the effective
degrees of freedom and the coefficients' own variance matrix do not
depend on it; the standard error of a hyperparameter and the Newton
direction of an outer search do depend on it.

The result is symmetric in `v` and `u` because mixed partial derivatives
are equal, and an implementation that pairs a direction with the wrong
parameter loses that symmetry.

The base method returns zeros, and that is exact. A design that does not
move with its coefficients has a second derivative that is identically
zero, which covers
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md),
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
and the five penalized constructors
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
when they have no method of their own.

For [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
the closed form is one order above the one that
[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md)
carries. With \\\theta_p = h_p(z_p)\\ and \\z_p = Z_p\beta\_{(p)}\\,
writing \\\tilde v_p = Z_p v\_{(p)}\\ and \\\tilde u_p = Z_p u\_{(p)}\\
for the directions carried onto each parameter's own scale,
\$\$\Big(\frac{\partial^2 X}{\partial\beta^2}\[v, u\]\Big)\[i, c_1\] =
Z\_{p_1}\[i, c_1\] \sum\_{p_2}\sum\_{p_3} r\_{p_1p_2p_3}(i)\\ \tilde
v\_{p_2}(i)\\\tilde u\_{p_3}(i),\$\$ \$\$r\_{p_1p_2p_3} =
f\_{p_1p_2p_3}h'\_{p_1}h'\_{p_2}h'\_{p_3} + \delta\_{p_1p_3}
f\_{p_1p_2}h''\_{p_1}h'\_{p_2} + \delta\_{p_2p_3}
f\_{p_1p_2}h'\_{p_1}h''\_{p_2} + \delta\_{p_1p_2}
f\_{p_1p_3}h''\_{p_1}h'\_{p_3} + \delta\_{p_1p_2p_3}
f\_{p_1}h'''\_{p_1}.\$\$ The five addends come from differentiating the
two addends of `term_block_deriv`'s \\q\_{p_1p_2}\\ once more: the first
from \\f\_{p_1p_2}\\, the second and third from \\h'\_{p_1}\\ and
\\h'\_{p_2}\\, the fourth from \\f\_{p_1}\\ inside the term that the
Kronecker delta carries, and the fifth from \\h''\_{p_1}\\. Exchanging
\\p_2\\ and \\p_3\\ sends the second addend to the fourth and leaves the
other three where they are, which is the symmetry in `v` and `u`.

This order needs no new derivation: \\f\_{p_1p_2p_3}\\ is the third
order of
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md),
served by the same routes as the lower orders, and \\h'''\\ is
[`linkfunctions7::d3linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d3linkinv.html),
exact for every shipped link and numerical for a user-defined one. The
cost is \\O(nP^3)\\, where \\P\\ is the number of the term's own
parameters, and the method is called once per pair of hyperparameters
rather than once per observation.

The result for a break-point term depends on its construction. If
`smoothed` is set to an
[`numericals7::abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.html),
the block is the true Jacobian and the closed forms are those of the
smoother, taken one order higher than
[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md)
reads them: with \\u = x - \psi\\, \\P = s''/2\\, \\T = s'''/2\\ and \\Q
= s''''/2\\, the change columns carry \\P\\ (for \\\gamma\\) and \\T\\
(for \\\delta\\) times the product of the two break-point directions in
their own entries, and the break-point column carries the two mixed
pieces and \\-(\gamma T + \delta Q)\\ times that product. The sharp
constructions return zeros: for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
the second derivative is genuinely zero away from the break-points,
because the truncated line's derivative in the position is an indicator
and the position column is linear in the change, while for
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
the block is a working linearization with a frozen weight, which is why
the first-order generics already return zeros there.

In the smoothed branch, where a break-point sits against its confinement
limit, the whole contribution is exactly zero, because every addend
carries a direction in the break-point and the confinement gates that
direction;
[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md)
is zero there for the same reason. Under
[`numericals7::smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.html),
which is exact outside \\\[-h, h\]\\, the answer is exactly zero on
every observation further than the width from a break-point. That
smoother is \\C^3\\, so its fourth derivative jumps at \\\pm h\\, and
the answer is exact at every observation that is not at distance \\h\\
from a break-point.

## See also

[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md),
[`term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.md),
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
dim(term_block_deriv2(b, v = c(1, 0), u = c(0, 1)))
#> [1] 20  2

# The result is symmetric in the two directions.
max(abs(term_block_deriv2(b, v = c(1, 0), u = c(0, 1)) -
        term_block_deriv2(b, v = c(0, 1), u = c(1, 0))))
#> [1] 0
```
