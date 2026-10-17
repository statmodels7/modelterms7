# The Hyperparameter Names That a Penalty Factory Produces

Builds the penalty at a probe size and reads its `params`, returning
`NULL` where it cannot be built.

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

A factory is a function of the coefficient count, and the names that it
produces do not depend on that count, so one coefficient is enough to
read them. The probe is wrapped because a factory may reject a size of
one, and the two outcomes are kept apart: `NULL` marks a penalty whose
names cannot be read here, and an empty vector marks a penalty with no
hyperparameter. Every shipped factory returns a penalty, so `NULL`
arises only from a factory written by a caller.
