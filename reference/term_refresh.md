# Refresh a Term at New Coefficients

Recomputes the quantities on which a term's design block depends, when
the coefficients move. For every ordinary term this is the identity,
because the block is a function of the data only. For a nonlinear term
the block is the Jacobian of its contribution, which depends on where
the parameters currently are.

## Usage

``` r
term_refresh(term, coef, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients at which to refresh, one per column of the term's
  block. For a nonlinear term these are the coefficients of its
  parameters' submodels, not the parameters themselves.

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
coefficient that it multiplies is the increment \\h\\. Refreshing at the
new \\\theta\\ and solving the linear problem again is the Gauss-Newton
iteration.
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
reports \\f(x; \theta)\\ itself, which is the other half that a step
needs.

A break-point term has the same shape with a different block: the
Jacobian for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
and for
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
the frozen-weight columns, from which the break-point is read off two
coefficients instead of being incremented.
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
distinguishes the two.

The fitting layer commits the refresh once per sweep and not at each
trial point. A discontinuous term's rescaling factor halves whenever the
break-point reverses direction, and is not allowed below a floor derived
from the conditioning of the working block, which the refresh records.
The factor depends on the path of the iteration: advancing it inside a
line search would anneal on every trial, and refreshing from the
specification each time would freeze the schedule at its starting value.

For a nonlinear term `coef` must have one entry per coefficient of the
term; a vector of another length is rejected.

## See also

[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
for the terms that implement it,
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
for the contribution that a step needs beside the block,
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
