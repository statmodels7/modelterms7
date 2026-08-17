# The Derivatives of a Nonlinear Term in Its Own Parameters

\\\partial^k f/\partial\theta^k\\ at the given coefficients, one named
component per index multiset, keyed as distributions7 keys its own
derivative surfaces.

## Usage

``` r
nl_fderiv(term, coef = NULL, order = 1L)
```

## Arguments

- term:

  A built
  [`NlTerm`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).

- coef:

  The coefficients, or `NULL` for the ones the term carries.

- order:

  1, 2, 3 or 4.

## Value

A named list of numeric vectors, one per component of that order.

## Details

The route is chosen per ORDER, highest first: a function supplied to
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md) is
used where there is one, then the symbolic route where the expression
can be differentiated, and otherwise one stencil applied to the highest
order that is analytic – which includes a supplied one. That is what
makes writing out a Hessian pay twice: the third and fourth orders are
then one difference away from an exact second rather than from the
function.

The components are in the term's OWN parameters, not in its
coefficients. The chain rule onto the coefficients – the links, and a
subformula's design – belongs to the term, which is the only thing that
knows them.

## See also

[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md)

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
names(nl_fderiv(b, order = 2))
#> [1] "a_a" "a_r" "r_r"
```
