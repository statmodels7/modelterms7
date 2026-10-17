# Where a Nonlinear Term's Coefficients Begin

Returns the coefficients at which
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
built the block: the starting value of each of the term's own parameters
carried through its link, from `start` where the caller gave one, and
coefficient zero for a parameter that `start` does not name. The block
is the Jacobian at those values, so it is a different block at any other
point. Where `target` is given (the response on the scale of the
predictor, which the fitting layer supplies), the parameters that
`start` does not fix are estimated from it by least squares. Zero
coefficients are a poor start for a nonlinear term, since the block is
then linearized at a point where the function was not meant to be
evaluated, and every quantity read at those coefficients, such as the
top of the path of a kinked penalty, depends on that point.

The search is a deterministic grid on each free parameter's own chart,
polished from the best grid points by
[`stats::optimize()`](https://rdrr.io/r/stats/optimize.html) for one
parameter and by Nelder-Mead for several, so that a start does not
depend on the random seed. The parameters in which the function is
jointly affine are solved by least squares at each point of the grid
instead of being searched over. They are found with
[`stats::D()`](https://rdrr.io/r/stats/deriv.html) by removing
repeatedly from the set of parameters every one whose derivative names a
member of the set (including the parameter itself); for a function given
as an opaque `fn` the set is empty. Where a parameter is developed over
a break-point term, each break-point is then started as described in
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md).

The box that each grid spans comes from the data: the quantiles, the
standard deviation and the range of the response and of every covariate
that the formula names, carried through the parameter's own link and
widened, because a parameter may be measured in the units of the
response or of a covariate and the box has to cover both.

It is a starting value and is allowed to be approximate; anything that
fails leaves the coefficients that
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
computed.

## Arguments

- term:

  A built
  [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).

- target:

  Optional numeric vector, one per observation: the response on the
  scale of the predictor. If it is absent, the coefficients are those
  that
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  computed.

- ...:

  Unused.

## Value

A numeric vector, one value per column of the block.
