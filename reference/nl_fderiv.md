# The Derivatives of a Nonlinear Term in Its Own Parameters

Returns \\\partial^k f/\partial\theta^k\\ at the given coefficients,
with one named component per index multiset, keyed as distributions7
keys its own derivative surfaces.

## Usage

``` r
nl_fderiv(term, coef = NULL, order = 1L)
```

## Arguments

- term:

  A built
  [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).
  Any other object is rejected with the error
  `"'term' must be a built nl() term."`, and an unbuilt term with
  `"the term is not built."`.

- coef:

  The coefficients at which to evaluate the derivatives, one per
  coefficient of the term, or `NULL` for the coefficients that the term
  currently carries, which
  [`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
  last committed. A vector of another length is rejected.

- order:

  1, 2, 3 or 4. Any other value is rejected.

## Value

A named list of numeric vectors of length `n`, one per component of that
order: `a`, `r` at order one, `a_a`, `a_r`, `r_r` at order two, and so
on, the parameter names joined by `"_"` in the term's own parameter
order.

## Details

The route is chosen for each order separately, in this priority: a
function supplied to
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
where there is one, then the symbolic route where the expression can be
differentiated, and otherwise one stencil applied to the highest
analytic order below the order wanted (a supplied order is included), or
to the value of \\f\\ if every lower order is computed numerically. A
supplied Hessian therefore serves the third and fourth orders, which are
then one difference away from an exact second order.

In each parameter the difference is a single stencil of the order
required, so a derivative of order four in one parameter is never built
from nested first differences. Stencils in different parameters are
composed.

The components are in the term's own parameters, not in its
coefficients. The chain rule onto the coefficients, through the links
and the design of a subformula, is applied by the term.

## See also

[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) for
the four functions that may supply these,
[`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md)
for the same quantities chained onto the coefficients.

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)

# Each order has one component per index multiset.
names(nl_fderiv(b, order = 1))
#> [1] "a" "r"
names(nl_fderiv(b, order = 2))
#> [1] "a_a" "a_r" "r_r"
names(nl_fderiv(b, order = 3))
#> [1] "a_a_a" "a_a_r" "a_r_r" "r_r_r"

# On a formula they are symbolic, so they are exact.
d1 <- nl_fderiv(b, order = 1)
e <- exp(-1.3 * dd$x)
c(a = max(abs(d1$a - e)), r = max(abs(d1$r + 2 * dd$x * e)))
#> a r 
#> 0 0 

# Order one is the Jacobian in the parameters, so with identity links
# and no subformula it is the design block itself.
max(abs(cbind(d1$a, d1$r) - term_matrix(b)))
#> [1] 0
```
