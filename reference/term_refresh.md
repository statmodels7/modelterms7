# Refresh a Term at New Coefficients

Recomputes whatever a term's design block depends on when the
coefficients move. For every ordinary term this is the identity, the
block being a function of the data alone. For a nonlinear term the block
is the Jacobian of its contribution, which depends on where the
parameters currently are.

## Usage

``` r
term_refresh(term, coef, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients to refresh at, as long as the term's block is wide.
  For a nonlinear term these are the coefficients of its parameters'
  submodels, not the parameters themselves.

- ...:

  Passed to methods.

## Value

A term of the same class with its block, and whatever else moves with
the coefficients, recomputed at `coef`. The base method returns `term`
unchanged.

## Details

A term whose contribution \\f(x; \theta)\\ is nonlinear in its
parameters is carried by the linearization

\$\$f(x; \theta + h) \approx f(x; \theta) + J(\theta) h, \qquad
J(\theta)\_{ij} = \frac{\partial f(x_i; \theta)}{\partial \theta_j},\$\$

so the design block at the current \\\theta\\ is \\J(\theta)\\ and the
coefficient it multiplies is the increment \\h\\. Refreshing at the new
\\\theta\\ and solving the linear problem again is the Gauss-Newton
iteration.
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
reports \\f(x; \theta)\\ itself, which is the other half a step needs.

A break-point term has the same shape with a different block: the
Jacobian for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
and for
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
the frozen-weight columns, from which the break-point is read off two
coefficients instead of being incremented.
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
is what tells the two apart.

The refresh is committed once per sweep, never per trial point. A
discontinuous term's rescaling factor halves whenever the break-point
reverses direction, which is a fact about the path: advancing it inside
a line search would anneal on every trial, and refreshing from the
specification each time would freeze the schedule at its starting value.

## See also

[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
for the terms that implement it,
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
for the contribution a step needs beside the block,
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
for which kind of block came back,
[`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
for the verdict on such a term.

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
built <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 1)), dd)
r1 <- term_refresh(built, c(2, 0.5))
max(abs(term_matrix(r1) - term_matrix(built))) > 0
#> [1] TRUE
```
