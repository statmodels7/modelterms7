# S7 Class for Nonlinear Parametric Terms

A subclass of
[`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a parametric function that is nonlinear in its own parameters. The
design block is the Jacobian of that function in the parameters, so the
term is linear in the sense the model layer needs while the function is
not; the block depends on where the parameters currently are and is
recomputed by
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md).

## Usage

``` r
NlTerm(
  label = character(0),
  hyper = list(),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  fn = NULL,
  nl_params = character(0),
  links = list(),
  subformulas = list(),
  deriv_mode = character(0),
  spec = list()
)
```

## Arguments

- label:

  A character string prefixed to the coefficient names when non-empty.

- hyper:

  The hyperparameters of the term's penalty that the caller HELD, as a
  named list. Empty, the default, means every one of them is estimated.
  See
  [`term_hyper`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md).

- X:

  The design block, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- coef_names:

  The coefficient names, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- blueprint:

  The information needed to reproduce the mapping on new data, filled by
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- penalty:

  The penalty on the term's coefficients, or `NULL`.

- fn:

  The function or formula defining the contribution.

- nl_params:

  The names of the nonlinear parameters.

- links:

  One link per parameter.

- subformulas:

  One optional formula per parameter.

- deriv_mode:

  How the derivatives are obtained.

- spec:

  The resolved construction settings.

## Value

An object of class `NlTerm`.

## See also

[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md)

## Examples

``` r
S7::S7_inherits(nl(~ theta1 * exp(theta2 * x)), NlTerm)
#> [1] TRUE
```
