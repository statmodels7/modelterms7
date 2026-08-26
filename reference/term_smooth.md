# Whether a Term's Penalized Objective Is Smooth

`TRUE` when the term's contribution to the penalized objective is
differentiable in the coefficients. A fitting layer reads it to split
the coefficient vector into the block a classical optimizer handles and
the block that needs a proximal step or a coordinate descent.

## Usage

``` r
term_smooth(term, ...)
```

## Arguments

- term:

  An object inheriting from
  [`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md).
  Build it first, or the answer is `TRUE` by default.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A single logical, never `NA`.

## The answer comes from the penalties

The term does not declare it. Every entry of
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
is asked for its kink set through
[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html),
at a probe value inside each hyperparameter's bounds, and the answer is
`FALSE` as soon as one entry reports a point. So an unpenalized term is
smooth, a ridge or a Gaussian prior is smooth, and lasso, SCAD, MCP and
the elastic net are not. A term cannot disagree with its own penalties.

The probe is any admissible value, the kink set being structural: the
midpoint of a bounded interval, one step inside a half-bounded one, zero
where the interval is the whole line.

## It answers for the whole term

A term carrying a penalty over part of its parameters and none over the
rest answers for the part. `nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g))`
is not smooth: the coefficients developing `a` sit at a kink, although
`r` is unpenalized.

## A specification is always smooth

The penalty is attached at
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
so
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
on an unbuilt term is empty and this answers `TRUE` whatever penalty the
term will carry. `term_smooth(lasso(~ x))` is `TRUE` and
`term_smooth(term_build(lasso(~ x), d))` is `FALSE`. Ask a built term.

## See also

[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the entries it runs over,
[`penalties7::penalty_kinks()`](https://statmodels7.github.io/penalties7/reference/penalty_kinks.html)
for the set it reads,
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md),
which uses the same split to count degrees of freedom.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(30), g = factor(rep(c("a", "b"), 15)))

# Unpenalized and quadratically penalized terms are smooth.
vapply(list(linpar(~ x), ridge(~ x), s(x, k = 5), random(~ 1 | g)),
       function(t) term_smooth(term_build(t, d)), logical(1))
#> [1] TRUE TRUE TRUE TRUE

# The four kinked penalties are not.
vapply(list(lasso(~ x), scad(~ x), mcp(~ x), enet(~ x)),
       function(t) term_smooth(term_build(t, d)), logical(1))
#> [1] FALSE FALSE FALSE FALSE

# A kink on part of a term's parameters makes the term non-smooth.
d$y <- 2 * exp(-1.3 * d$x)
nb <- term_build(nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g),
                    start = list(r = 1.3)), d)
term_smooth(nb)
#> [1] FALSE
c(npar = term_npar(nb), penalized = length(term_penalties(nb)[[1]]$index))
#>      npar penalized 
#>         3         2 

# Unbuilt, there is no penalty to read yet.
c(spec = term_smooth(lasso(~ x)), built = term_smooth(term_build(lasso(~ x), d)))
#>  spec built 
#>  TRUE FALSE 
```
