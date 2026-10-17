# Defining a model term

In many fitting functions a term in a model formula is a string on which
the function switches. Writing a new term then means editing that
function, and a user cannot write one at all. Here a term is an object.
It builds its own design block, carries its own penalty, reports its
number of parameters and reapplies itself to new data. A formula
interpreter recognizes a term by the value that its call *returns*, so a
term defined outside the package works in a formula with no edit
anywhere.

This vignette writes a term that the package does not ship: a seasonal
term, a Fourier design of $`2k`$ columns over a known period, penalized
so that the higher harmonics shrink first.

## The minimal additive term

A new additive term is a subclass of `additive_term` together with a
method for
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md).
State that belongs to the term is declared as properties of its own
subclass; the shipped `SmoothTerm` carries `vars`, `by`, `spec` and
`sparse` in this way.

``` r

SeasonTerm <- S7::new_class(
  "SeasonTerm", parent = additive_term,
  properties = list(expr = S7::class_any, spec = S7::class_list)
)
```

The constructor captures the covariate unevaluated and records the
settings that the build needs. The constructor computes nothing, because
a specification is not a block and the data are not yet in scope.

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
`blueprint` and `penalty`. The blueprint stores what
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
needs later, so that prediction reapplies recorded quantities and does
not recompute them from the new rows. Storing the blueprint at build
time is the second requirement on a new term.

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

The penalty is the quadratic form
$`\tfrac{\lambda}{2}\sum_j j^2(a_j^2 + b_j^2)`$, plus the normalizing
constant of the corresponding Gaussian prior, where $`a_j`$ and $`b_j`$
are the cosine and sine coefficients of harmonic $`j`$ and $`\lambda`$
is the smoothing parameter. Harmonic $`j`$ oscillates $`j`$ times as
fast, so the weight $`j^2`$ penalizes the higher harmonics more
strongly, and $`\lambda`$ controls the smoothness of the fit as well as
its size.

The constructor and the
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
method are the whole of the compulsory part. The term can now be built
on data:

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

A least-squares fit on the block recovers the coefficients that
generated the data:

``` r

cf <- lm.fit(cbind(1, term_matrix(tm)), d$y)$coefficients
names(cf) <- c("(Intercept)", term_coef_names(tm))
round(cf, 3)
#> (Intercept)        cos1        cos2        cos3        sin1        sin2 
#>       0.011      -0.019      -0.022       0.369       2.049      -0.001 
#>        sin3 
#>       0.001
```

The truth is `sin1 = 2`, `cos3 = 0.4` and zero for the other
coefficients.

## Prediction at new data

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
has no base method, so a term that is used for prediction or for
cross-validation defines one. The method reapplies the blueprint:

``` r

S7::method(term_predict, SeasonTerm) <- function(term, newdata, ...) {
  x <- as.numeric(eval(term@expr, newdata, baseenv()))
  bp <- term@blueprint
  .season_X(x, bp$period, bp$k)
}
```

The period comes from the blueprint and never from the new data. A term
that recomputed it would be a different term on every subset. The
identity that detects this mistake is that prediction on the rows of the
fit must return the block of the fit:

``` r

rows <- c(3, 40, 111)
all.equal(term_predict(tm, d[rows, , drop = FALSE]),
          term_matrix(tm)[rows, , drop = FALSE], check.attributes = FALSE)
#> [1] TRUE
```

For a seasonal term the period is a constant that the caller supplied,
so the hazard is small. For a spline the recorded quantities are the
knots, and for a factor they are the levels. Each is learned from the
fitting data, and each would silently change if it were learned again at
prediction.

The [`baseenv()`](https://rdrr.io/r/base/environment.html) enclosure is
deliberate. The expression is evaluated in the new data with only the
base namespace behind it, so a variable of the calling environment is
not found. Under any other enclosure such a variable would be picked up
from the calling environment and would quietly supply the fitting data’s
values.

## Validating the term

[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
builds the term and runs six checks on the result:

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

The last check builds the term on all the rows, calls
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
on a subset of them and compares the result with the corresponding rows
of the block. A term that rebuilds its basis or relearns its levels at
prediction fails this check. When the data carry a factor, the subset
omits every row of the factor’s last level; otherwise it is the first
half of the rows, as here. The subset is passed through
[`droplevels()`](https://rdrr.io/r/base/droplevels.html) first, because
a plain row subset of a data frame keeps every factor level, including
levels that do not occur in any row, and a term that rebuilt its levels
would look correct on it.

## The penalty and its rank

The penalty travels with the term, so the effective degrees of freedom
follow from the term alone:

``` r

H <- crossprod(term_matrix(tm))
c(`lambda = 1` = edf(tm, cf[-1], H, theta = list(lambda = 1)),
  `lambda = 1e6` = edf(tm, cf[-1], H, theta = list(lambda = 1e6)))
#>   lambda = 1 lambda = 1e6 
#> 5.7381084609 0.0002722007
```

The count runs from 6 down to 0, because the penalty above has full
rank: every coefficient is penalized, so a large smoothing parameter
removes the seasonal effect entirely. A shipped
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) with
its default smoother has a penalty that is rank deficient by exactly
one, so that the effective count runs from the number of coefficients
down to one, and the limit is a straight line. The cyclic, Fourier and
first-order smoothers have no free line, and their count goes to zero.
Which of the two behaviors a term has depends on the null space of its
penalty, and the choice of that null space is a modeling decision.

[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
reports the entries that a term declares. Each entry names a subset of
the term’s own parameters and the penalty over them, together with the
hyperparameters that are held and the settings of the path over the
others. One entry that covers every parameter is the ordinary case:

``` r

length(term_penalties(tm))
#> [1] 1
names(term_penalties(tm)[[1]])
#> [1] "name"      "index"     "penalty"   "fixed"     "n_values"  "values"   
#> [7] "min_ratio" "search"    "ids"
```

A term with several penalties declares several entries. A
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
whose developments carry a penalty on each kind of coefficient is an
example: one entry covers the changes of slope and another the changes
of level. A slope change and a jump are not comparable quantities, and
one hyperparameter over their union would price them against each other.
A
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
written without developments declares no entry.

``` r

jj <- term_build(jseg(t, gamma ~ 0 + lasso(~ 1), delta ~ 0 + ridge(~ 1)), d)
vapply(term_penalties(jj), function(e) e$name, character(1))
#> [1] "gamma1::lasso(~1)" "delta1::ridge(~1)"
length(term_penalties(term_build(jseg(t), d)))
#> [1] 0
```

## In a formula

[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
recognizes a term by evaluating its call and testing whether the value
inherits `model_term`. The term names are never listed, so the new term
is recognized with no edit:

``` r

d$x <- rnorm(nrow(d))
f <- interpret_formula(y ~ x + season(t, period = 1), data = d)
vapply(f$terms, function(z) class(z)[[1]], character(1))
#>                    linpar     season(t, period = 1) 
#> "modelterms7::LinparTerm"              "SeasonTerm"
```

A bare covariate collapses into one parametric block. A call whose value
does not inherit `model_term` stays a covariate when a model matrix can
hold the value, so `log(x)` is arithmetic. A call whose evaluation fails
signals an error that names the call.

Recognition by evaluation also produces a third case, a masked name.
`mgcv` exports
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md), and
`segmented` exports
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
so a user who has attached either package and writes this package’s
formula gets the function of the other package. The interpreter rejects
the result at the call itself. The error names the term, the class of
the value and the package that supplied the function, and it suggests
`modelterms7::`, so the value does not reach
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) and fail
there.

## Inherited methods and structural terms

[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md),
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md),
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
[`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md),
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md),
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md),
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) are inherited
from the base classes and read the properties that
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
filled.

Two further generics matter only for a term whose block moves with its
coefficients.
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
recomputes the block at new coefficients.
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
reports whether that block is a Jacobian, so that the fitting step is
Gauss-Newton, or a frozen working linearization, so that the step is a
fixed-point iteration. For
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and every smoothed break-point term, the answer is `TRUE`. For a sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
the answer is `FALSE`, because its break-point is read off a product of
coefficients and is not estimated as a parameter. A term with a fixed
design, like the one above, inherits a
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
that returns the term unchanged and a
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
that returns `TRUE`, and it needs neither method.

The other branch of the class tree is `structural_term`, whose terms
rewrite the likelihood instead of adding columns:
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
drives a distribution parameter by a score-driven recursion, and
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
mixes the likelihood over a latent Markov chain. A structural term has a
different structure. It supplies
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
and either
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
or
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md),
and
[`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
carries the derivative of a filter.
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
fits at most one structural term in a formula.

## Summary

- Minimum to define an additive term: a subclass of `additive_term` and
  a
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  method that fills `X`, `coef_names`, `blueprint` and `penalty`. Any
  state of its own is a property of the subclass.
- The constructor captures the covariate unevaluated and computes
  nothing.
- [`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  has no base method. Write one, reapply the blueprint, and evaluate in
  the new data with
  [`baseenv()`](https://rdrr.io/r/base/environment.html) as the
  enclosure.
- [`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
  runs six checks, and the last one runs on a
  [`droplevels()`](https://rdrr.io/r/base/droplevels.html) subset, which
  exposes a rebuild.
- The null space of the penalty determines what a large smoothing
  parameter leaves behind, and
  [`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
  follows from the term alone.
- [`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
  recognizes a term by evaluating its call, so a term defined outside
  the package needs no registration.
- A block that moves with its coefficients adds
  [`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
  and
  [`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md);
  a term that rewrites the likelihood belongs to the other branch.
