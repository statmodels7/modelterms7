# Defining a model term

A term in a model formula is usually a string a fitting function
switches on. Writing a new one then means editing that function, and a
user cannot write one at all. Here a term is an object: it builds its
own design block, carries its own penalty, says how many parameters it
spends, and reapplies itself to new data. A formula interpreter
recognizes it by what its call *returns*, so a term defined outside the
package drops into a formula with no edit anywhere.

This vignette writes one the package does not ship: a seasonal term, a
Fourier design of $`2k`$ columns over a known period, penalized so that
the higher harmonics shrink first.

## The minimum

A new additive term is a subclass of `additive_term` and a method for
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md).
A term declares any state of its own as **properties on its own
subclass**; the shipped `SmoothTerm` carries `vars`, `by` and `spec`
that way.

``` r

SeasonTerm <- S7::new_class(
  "SeasonTerm", parent = additive_term,
  properties = list(expr = S7::class_any, spec = S7::class_list)
)
```

The constructor captures the covariate **unevaluated** and records what
the build will need. Nothing is computed here: a specification is not a
block, and the data are not in scope yet.

``` r

season <- function(x, period, k = 3, lambda = NULL) {
  xe <- substitute(x)
  SeasonTerm(
    label = sprintf("season(%s)", deparse(xe)),
    expr = xe,
    spec = list(period = period, k = k),
    hyper = if (is.null(lambda)) list() else list(lambda = lambda),
    X = NULL, coef_names = character(0),
    blueprint = list(), penalty = NULL
  )
}
```

[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
evaluates the expression in the data and fills `X`, `coef_names`,
`blueprint` and `penalty`. The **blueprint** is what
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
will need later, and putting it there rather than recomputing it is the
whole of the contract’s second half.

``` r

.season_X <- function(x, period, k) {
  j <- seq_len(k)
  w <- outer(2 * pi * x / period, j)
  X <- cbind(cos(w), sin(w))
  colnames(X) <- c(sprintf("cos%d", j), sprintf("sin%d", j))
  X
}

S7::method(term_build, SeasonTerm) <- function(term, data, ...) {
  x <- as.numeric(eval(term@expr, data, baseenv()))
  sp <- term@spec
  S7::set_props(
    term,
    X = .season_X(x, sp$period, sp$k),
    coef_names = c(sprintf("cos%d", seq_len(sp$k)),
                   sprintf("sin%d", seq_len(sp$k))),
    blueprint = sp,
    penalty = penalties7::quadratic_penalty(diag(rep(seq_len(sp$k)^2, 2)))
  )
}
```

The penalty is $`\sum_j j^2(a_j^2 + b_j^2)`$, the discrete analogue of a
roughness penalty: harmonic $`j`$ oscillates $`j`$ times as fast, so
weighting by $`j^2`$ makes the smoothing parameter buy smoothness, not
shrinkage alone.

That is the whole of the compulsory part.

``` r

d <- data.frame(t = seq(0, 4, length.out = 200))
d$y <- 2 * sin(2 * pi * d$t) + 0.4 * cos(6 * pi * d$t) + rnorm(200, sd = 0.3)

tm <- term_build(season(t, period = 1, k = 3), d)
c(built = term_is_built(tm), npar = term_npar(tm))
#> built  npar 
#>     1     6
term_coef_names(tm)
#> [1] "cos1" "cos2" "cos3" "sin1" "sin2" "sin3"
dim(term_matrix(tm))
#> [1] 200   6
```

and it recovers what it was given:

``` r

cf <- lm.fit(cbind(1, term_matrix(tm)), d$y)$coefficients
names(cf) <- c("(Intercept)", term_coef_names(tm))
round(cf, 3)
#> (Intercept)        cos1        cos2        cos3        sin1        sin2 
#>       0.011      -0.019      -0.022       0.369       2.049      -0.001 
#>        sin3 
#>       0.001
```

against a truth of `sin1 = 2`, `cos3 = 0.4` and zero elsewhere.

## Prediction reapplies; it does not rebuild

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
has **no base method**, so a term that will be predicted from or
cross-validated writes one, and what it must do is reapply the
blueprint:

``` r

S7::method(term_predict, SeasonTerm) <- function(term, newdata, ...) {
  x <- as.numeric(eval(term@expr, newdata, baseenv()))
  bp <- term@blueprint
  .season_X(x, bp$period, bp$k)
}
```

The period comes from the blueprint, never from the new data. A term
that recomputed it would be a different term on every subset, and the
identity that catches the mistake is that predicting on rows the model
was fitted to must return the block it was fitted with:

``` r

rows <- c(3, 40, 111)
all.equal(term_predict(tm, d[rows, , drop = FALSE]),
          term_matrix(tm)[rows, , drop = FALSE], check.attributes = FALSE)
#> [1] TRUE
```

For a seasonal term the period is a constant the caller gave, so the
hazard is mild. For a spline it is the knots, for a factor the levels,
and for a standardized penalized block the column spreads: each is
learned from the fitting data and each would silently move if it were
relearned.

The [`baseenv()`](https://rdrr.io/r/base/environment.html) enclosure is
deliberate. It means the expression is evaluated in the new data and
nowhere else, so a name that is not a column is not found, where it
would otherwise be picked up from the calling environment and quietly
supply the fitting data’s values.

## Validating it

[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
builds the term, then asks six questions of the result:

``` r

print(check_term(season(t, period = 1, k = 3), d, verbose = FALSE))
#>       check status          info
#> 1     build     OK 200 x 6 block
#> 2     names     OK              
#> 3      npar     OK              
#> 4    smooth     OK        smooth
#> 5 reproduce     OK              
#> 6    subset     OK      100 rows
```

The last is the one worth understanding. It builds the term on a
**subset** of the rows and compares against predicting there, which is
how a rebuilt basis or a relearned set of levels is exposed. The subset
is passed through
[`droplevels()`](https://rdrr.io/r/base/droplevels.html) first: a plain
row subset of a data frame keeps every factor level, including the ones
no row uses, so a term that rebuilt its levels would look correct on it.

## The penalty, and what its rank decides

The penalty travels with the term, so the effective degrees of freedom
follow from the term alone:

``` r

H <- crossprod(term_matrix(tm))
c(`lambda = 1` = edf(tm, cf[-1], H, theta = list(lambda = 1)),
  `lambda = 1e6` = edf(tm, cf[-1], H, theta = list(lambda = 1e6)))
#>   lambda = 1 lambda = 1e6 
#> 5.7381084609 0.0002722007
```

The count runs from 6 down to 0, because the penalty above has **full
rank**: nothing is left unpenalized, so a large smoothing parameter
removes the seasonal effect entirely. A shipped
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) is
deliberately different, its penalty being rank deficient by exactly one
so that the effective count runs from $`k`$ down to **one** and the
limit is a straight line, not nothing. Which of the two a term wants is
a modeling decision, and it is made by choosing the penalty’s null
space.

[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
reports the entries a term declares, each naming a subset of its own
parameters and the penalty over them. One entry covering everything is
the ordinary case:

``` r

length(term_penalties(tm))
#> [1] 1
names(term_penalties(tm)[[1]])
#> [1] "name"      "index"     "penalty"   "fixed"     "n_values"  "values"   
#> [7] "min_ratio" "search"
```

A term with several may declare several:
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
declares two, one on its slope changes and one on its jumps, because a
slope change and a jump are not comparable quantities and one
hyperparameter over their union would price them against each other.

## In a formula

[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
recognizes a term by **evaluating** its call and asking whether the
value inherits `model_term`. Nothing lists the term names, so the new
term is recognized with no edit:

``` r

d$x <- rnorm(nrow(d))
f <- interpret_formula(y ~ x + season(t, period = 1), data = d)
vapply(f$terms, function(z) class(z)[[1]], character(1))
#>                    linpar     season(t, period = 1) 
#> "modelterms7::LinparTerm"              "SeasonTerm"
```

A bare covariate collapses into one parametric block, and a call whose
value does not inherit `model_term` stays a covariate, so `log(x)` is
arithmetic.

⚠️ Recognition by evaluation has one consequence worth knowing: a
**masked name** is a third case. `mgcv` exports
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md),
`segmented` exports
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
so a user with either attached writes this package’s formula and gets
the other package’s function. The interpreter rejects that where it
happens, naming the class returned and the package that supplied the
function, instead of letting the value reach
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) and fail
there.

## What comes free, and the second branch

[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md),
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md),
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md),
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md),
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md),
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) all come from
the base classes and read the properties
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
filled.

Two further methods matter only for a term whose block **moves with its
coefficients**.
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
recomputes the block at new coefficients and
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
says whether that block is a Jacobian (so the fitting step is
Gauss-Newton) or a frozen working linearization (so it is a fixed-point
iteration).
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
answer `TRUE`;
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
answer `FALSE`, their break-point being read off a product of
coefficients, never estimated as one. A term with a fixed design, like
the one above, inherits the identity and needs neither.

The other branch of the class tree is `structural_term`, whose terms
rewrite the likelihood instead of adding columns:
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
drives a distribution parameter by a score-driven recursion, and
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
mixes the likelihood over a latent Markov chain. They are a different
contract, made of
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md),
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md),
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
and
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md),
and a formula may carry at most one of them.

## Summary

- **Minimum to define an additive term:** a subclass of `additive_term`
  and a
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  method that fills `X`, `coef_names`, `blueprint` and `penalty`. Any
  state of its own is a property on the subclass.
- The constructor captures the covariate unevaluated and computes
  nothing.
- [`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  has no base method. Write one, reapply the blueprint, and evaluate in
  the new data with
  [`baseenv()`](https://rdrr.io/r/base/environment.html) as the
  enclosure.
- [`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
  asks six questions, the last of them on a
  [`droplevels()`](https://rdrr.io/r/base/droplevels.html) subset, the
  one that exposes a rebuild.
- The penalty’s null space decides what a large smoothing parameter
  leaves behind, and
  [`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
  follows from the term alone.
- [`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
  recognizes a term by evaluating its call, so a term defined outside
  the package needs no registration.
- A block that moves with its coefficients adds
  [`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
  and
  [`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md);
  a term that rewrites the likelihood is the other branch entirely.
