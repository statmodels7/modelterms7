# Filter Several Score-Driven Terms as One Recursion

Checks that a set of built
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
terms can be filtered together and returns them as a
[`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md)
object, with the cross loadings that each term's `p_cross` declares
added to its parameters.

## Usage

``` r
gas_joint(terms)
```

## Arguments

- terms:

  A named list of built
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md)
  objects, named by the distribution parameter in whose equation each
  one sits.

## Value

A
[`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md)
object. Its `terms` are the input terms with `cross_with` filled, so
their
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
include the cross loadings.

## Details

The terms must walk the same series: their `by` and `time` must give the
same groups and the same order, because each score is read at the
predictors of every equation at the same time. A term may develop its
level, its loadings or its persistence with a subformula; its cross
loadings stay scalar.

The terms must share one `scaling`, because the driving vector is scaled
by a power of the expected information of all the filtered parameters
together. With more than one term the power is 0, 1/2 or 1: for these
the derivatives of the matrix power are exact and short (the identity,
the symmetric square root through a Sylvester equation per order, and
the inverse).

The `p_cross` of a term names the parameters whose scores its level
loads besides its own, with the number of lags of each. A parameter that
it names with a positive count must be one of the joined terms, and a
term cannot name its own parameter, whose score it loads through its own
loadings. The cross loadings of a term follow the parameters in the
order of the joined terms, lag by lag.

## See also

[`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md),
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = rnorm(40))
tm <- list(mu = term_build(gas(p = 1, q = 1, time = t,
                               p_cross = c(sigma = 2)), dd),
           sigma = term_build(gas(p = 1, q = 1, time = t), dd),
           nu = term_build(gas(p = 1, q = 1, time = t), dd))
jt <- gas_joint(tm)
lapply(jt@terms, term_params)
#> $mu
#> [1] "omega"        "kappa1"       "pacf1"        "kappa1_sigma" "kappa2_sigma"
#> 
#> $sigma
#> [1] "omega"  "kappa1" "pacf1" 
#> 
#> $nu
#> [1] "omega"  "kappa1" "pacf1" 
#> 
```
