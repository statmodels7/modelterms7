# The Second Derivative of a Design Block in the Coefficients

\\\sum\_{q,r} v_q u_r\\\partial^2
X\_{ij}/\partial\beta_q\partial\beta_r\\, one entry per observation and
column of the block: the block's second derivative contracted in two
directions the caller supplies.

## Usage

``` r
term_block_deriv2(term, coef = NULL, v, u, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients, or `NULL` for the ones it carries.

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
and it is contracted, never returned as an array, for the reason
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
is contracted in the structural branch: the full object has \\m\\
coefficient indices twice over beside the \\n\\ rows and \\m\\ columns
of the block, and only its contraction in the two directions the
penalized mode moves is ever read.

The quantity enters the Hessian of a marginal criterion and nothing
else. It is absent from the criterion's value, from its gradient and
from the fit, so the coefficients, the log-likelihood, the effective
degrees of freedom and the coefficients' own variance matrix do not
depend on it; what does depend on it is the standard error of a
hyperparameter and the Newton direction of an outer search.

The result is symmetric in `v` and `u`, mixed partial derivatives being
equal, and an implementation that pairs a direction with the wrong
parameter loses that symmetry.

The base method returns zeros, and that is exact. A design that does not
move with its coefficients has a second derivative that is identically
zero, which covers
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md),
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
and the five penalized constructors without a method of their own.

For [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
the closed form is one order above the one
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
\\h'\_{p_2}\\, the fourth from \\f\_{p_1}\\ inside the term the
Kronecker delta carries, and the fifth from \\h''\_{p_1}\\. Exchanging
\\p_2\\ and \\p_3\\ sends the second addend to the fourth and leaves the
other three where they are, which is the symmetry in `v` and `u`.

Nothing new is derived: \\f\_{p_1p_2p_3}\\ is the third order of
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md),
served by the same four-way machinery as the lower orders, and \\h'''\\
is
[`linkfunctions7::d3linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d3linkinv.html),
exact for every shipped link and numerical for a user-defined one. The
cost is \\O(nP^3)\\ in the term's own parameters, of which there are two
to four in practice, and the call is made once per pair of
hyperparameters rather than once per observation.

A break-point term answers according to its construction. With
`smoothed` an
[`penalties7::abs_smoother()`](https://statmodels7.github.io/penalties7/reference/abs_smoother.html)
the block is the true Jacobian and the closed forms are the smoother's
own one order further up than
[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md)
reads them: with \\u = x - \psi\\, \\P = s''/2\\, \\T = s'''/2\\ and \\Q
= s''''/2\\, the change columns contribute \\P\\ and \\T\\ twice in the
break-point and the break-point column contributes the two mixed pieces
and \\-(\gamma T + \delta Q)\\. The sharp constructions answer zeros:
for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
the second derivative is genuinely zero away from the break-points, the
truncated line's derivative in the position being an indicator and the
position column being linear in the change, while for
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
the block is a working linearization with a frozen weight, which is why
the first-order generics already answer zeros there.

Two properties of the smoothed branch are worth stating because they are
exact. Where a break-point sits against its confinement limit the whole
contribution is zero, every addend carrying a direction in the
break-point, which the first derivative does not, the position column
moving with the change whatever the position does. And under
[`penalties7::smooth_quintic()`](https://statmodels7.github.io/penalties7/reference/smooth_quintic.html),
which is exact outside \\\[-h, h\]\\, the answer is zero on every
observation further than the width from a break-point; that smoother is
\\C^3\\, so its fourth derivative jumps at \\\pm h\\ and the answer is
exact away from those two points, leaving it exact everywhere else.

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

# symmetric in the two directions
max(abs(term_block_deriv2(b, v = c(1, 0), u = c(0, 1)) -
        term_block_deriv2(b, v = c(0, 1), u = c(1, 0))))
#> [1] 0
```
