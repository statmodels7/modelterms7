# S7 Class for the Jeffreys Penalty

`JeffreysPenalty` is the subclass of
[`penalties7::penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.html)
that a
[`jfp()`](https://statmodels7.github.io/modelterms7/reference/jfp.md)
term attaches to its block. It carries the power of the determinant and
the setting of the constant direction, and it has no hyperparameters.

## Usage

``` r
JeffreysPenalty(
  penalty_name = character(0),
  map = NULL,
  n_coef = integer(0),
  params = character(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  power = integer(0),
  intercept = NULL
)

jeffreys_penalty(n_coef, power = 0.5, intercept = NULL)
```

## Arguments

- penalty_name, map, params, params_bounds, link_params, params_smooth:

  Properties inherited from
  [`penalties7::penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.html).
  The penalty has no map and no hyperparameters, so `jeffreys_penalty()`
  sets the name `"jeffreys"`, `map = NULL`, an empty `params`, empty
  lists of bounds and links, and an empty logical `params_smooth`.

- n_coef:

  The number of coefficients in the block.

- power:

  The power \\a\\ of the determinant, a single positive number.

- intercept:

  `NULL`, `TRUE` or `FALSE`, the setting of the constant direction (see
  [`jfp()`](https://statmodels7.github.io/modelterms7/reference/jfp.md)).

## Value

`jeffreys_penalty()` returns an S7 object of class `JeffreysPenalty`,
inheriting from
[`penalties7::penalty()`](https://statmodels7.github.io/penalties7/reference/penalty.html).

## Details

The penalty is \\-a\log\|\tilde I(\beta)\|\\, where \\\tilde I\\ is the
expected Fisher information of the block's coefficients, augmented by
the constant direction of the formula level where
[`jfp()`](https://statmodels7.github.io/modelterms7/reference/jfp.md)
includes it. That matrix depends on the family and on every predictor of
the model, so the penalty cannot be evaluated from the block's
coefficients alone, and the methods of
[`penalties7::penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.html),
[`penalties7::penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.html),
[`penalties7::penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.html)
and of the other generics that read the coefficients signal an error.
The fitting layer evaluates it on the model.

The methods that do not need the model return their answer:
[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
is empty,
[`penalties7::is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.html)
and
[`penalties7::is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.html)
are `FALSE`.

## See also

[`jfp()`](https://statmodels7.github.io/modelterms7/reference/jfp.md),
the term that attaches it.

## Examples

``` r
p <- jeffreys_penalty(3)
p@power
#> [1] 0.5
p@params
#> character(0)
penalties7::penalty_kinks(p, list())
#> numeric(0)
penalties7::is_proper(p)
#> [1] FALSE
```
