# The Hyperparameter Names a Penalty Factory Produces

Builds the penalty at a probe size and reads its `params`, returning an
empty vector where it cannot be built.

## Usage

``` r
.factory_params(factory)
```

## Arguments

- factory:

  A function of the coefficient count returning a penalty, or a penalty.

## Value

A character vector, empty where the penalty carries no hyperparameter,
or `NULL` where it could not be built.

## Details

A factory is a function of the coefficient count, and the names it
produces do not depend on that count, so one coefficient is enough to
read them. The probe is wrapped because a factory may reject a size of
one, and the two outcomes are kept apart: `NULL` says nothing is known
here, an empty vector says the penalty carries no hyperparameter. Every
shipped factory answers, so `NULL` is reachable only from one a caller
wrote.
