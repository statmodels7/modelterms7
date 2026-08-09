# Refresh a Term at New Coefficients

Recomputes whatever a term's design block depends on when the
coefficients move. For every ordinary term this is the identity, the
block being a function of the data alone; for a nonlinear term the block
is the Jacobian of its contribution, which is a function of where the
parameters currently are.

## Usage

``` r
term_refresh(term, coef, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The current coefficients of the term's block.

- ...:

  Passed to methods.

## Value

A built term, refreshed.

## Details

A term whose contribution \\f(x; \theta)\\ is nonlinear in its
parameters is carried by the linearization

\$\$f(x; \theta + h) \approx f(x; \theta) + J(\theta) h, \qquad
J(\theta)\_{ij} = \frac{\partial f(x_i; \theta)}{\partial \theta_j},\$\$

so the design block at the current \\\theta\\ is \\J(\theta)\\ and the
coefficient it multiplies is the increment \\h\\. Refreshing at the new
\\\theta\\ and solving the linear problem again is the Gauss-Newton
iteration;
[`term_value`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
reports \\f(x; \theta)\\ itself, which is the other half a step needs.
For a break-point term the same shape holds with a different block: the
Jacobian for
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md), and
the frozen-weight columns of
[`jump`](https://statmodels7.github.io/modelterms7/reference/seg.md),
from which the break-point is read rather than incremented.

## See also

[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md),
[`term_value`](https://statmodels7.github.io/modelterms7/reference/term_value.md)

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
built <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 1)), dd)
r1 <- term_refresh(built, c(2, 0.5))
max(abs(term_matrix(r1) - term_matrix(built))) > 0
#> [1] TRUE
```
