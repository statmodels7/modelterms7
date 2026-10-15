# The Smooth on the Left Side of a Bar

Reads the left side of a bar formula and returns the smooth written
there, with whether the groups also carry a random level, or `NULL`
where the left side holds no
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) or
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)
call.

## Usage

``` r
.random_smooth_lhs(within, env)
```

## Arguments

- within:

  The left side of the bar, a language object.

- env:

  The formula's environment.

## Value

`NULL`, or a list with `term`, the unbuilt
[`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md),
and `intercept`, a single logical.

## Details

The left side is split on `+`. A smooth random effect holds exactly one
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) or
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)
call, optionally beside `1`, `0` or `-1`, which follow the intercept
convention of a model formula. Any other covariate beside the smooth
signals an error naming the spelling that writes it as a term of its
own.

The call is evaluated with this package's
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) in
place of whatever the formula's environment binds to those names, so
that a session with mgcv attached builds the same term.

The smooth may not carry a `by` variable, a `sparse` setting, held
hyperparameters or an `id`. The grouping is the variable after the bar,
and the hyperparameters belong to the random effect, where `hyper` and
`id` of
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
name them.
