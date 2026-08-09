# modelterms7 <img src="man/figures/logo.png" align="right" height="139" alt="" />

Model terms as S7 objects. A term is the recipe a formula names: the
mapping from a data frame to a design-matrix block, the penalty attached
to the block's coefficients, and the metadata a fit reads -- coefficient
names, hyperparameters with their links, whether the penalized objective
is differentiable, and how effective degrees of freedom are counted.

The formula interpreter recognizes a term by what its call evaluates to,
not by its name, so a term class defined outside the package works in a
formula without registration. `linpar()` is the unpenalized parametric
block; `ridge()`, `lasso()`, `scad()` and `mcp()` attach the
corresponding [penalties7](https://statmodels7.github.io/penalties7/)
object to their block; `random()` builds grouped intercepts and slopes
with the effect distribution as the penalty, correlated or not through a
[parameters7](https://statmodels7.github.io/parameters7/) structure
replicated across groups; `s()` and `te()` are the penalized smooths of
one and of several covariates, built on
[basis7](https://statmodels7.github.io/basis7/); `nl()` carries a
contribution nonlinear in its own parameters, its block the Jacobian
refreshed as they move; `seg()`, `jump()` and `jseg()` estimate the
break-points at which an effect changes slope, level, or both; `gas()`
and `regime()` rewrite the likelihood rather than adding to it, with
score-driven dynamics and a latent Markov chain; `cens()` marks a
censored response on the left side of the formula. Part of the
[statmodels7](https://statmodels7.github.io) toolkit.

## Installation

``` r
# install.packages("pak")
pak::pak("statmodels7/modelterms7")
```

## A formula, interpreted

``` r
library(modelterms7)

dd <- data.frame(y = rnorm(50), x1 = rnorm(50), x2 = rnorm(50),
                 g = factor(sample(letters[1:5], 50, TRUE)))
dd$L <- matrix(rnorm(50 * 3), 50, 3)

out <- interpret_formula(y ~ x1 + log(abs(x2)) + lasso(L) +
                           random(~ 1 | g), dd)
names(out$terms)

built <- lapply(out$terms, term_build, data = dd)
sapply(built, term_npar)
sapply(built, term_smooth)
```

Every built term carries its design block, its coefficient names, and a
blueprint that reproduces the mapping on new data; the penalized terms
carry their [penalties7](https://statmodels7.github.io/penalties7/)
object, whose hyperparameters, derivatives and kink set the model layer
reads. `check_term()` validates a term against a data frame, and `edf()`
computes effective degrees of freedom from the pieces a fit supplies.
