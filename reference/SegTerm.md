# S7 Class for Break-Point Terms

The subclass of
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
for a covariate whose effect changes at estimated break-points: a change
of slope
([`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)),
a change of level
([`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)),
or both at the same points
([`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)).
Its design block is the **working** block of the iteration that
estimates the break-points, and
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
recomputes it as they move.

## Usage

``` r
SegTerm(
  label = character(0),
  hyper = list(),
  grid = list(),
  values = list(),
  min_ratio = numeric(0),
  search = character(0),
  X = NULL,
  coef_names = character(0),
  blueprint = list(),
  penalty = NULL,
  kind = character(0),
  var = NULL,
  npsi = integer(0),
  linear = logical(0),
  subformulas = list(),
  spec = list()
)
```

## Arguments

- label:

  A character string prefixed to the term's coefficient names when
  non-empty, and used as the title of
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and the tag
  of [`print()`](https://rdrr.io/r/base/print.html). `character(0)` and
  `""` both mean no label.

- hyper:

  The hyperparameters of the term's penalty the caller **held**, as a
  named list keyed by the penalty's names. Empty, the default, means
  every one of them is estimated. See
  [`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md).

- grid:

  How many values a path visits for each of the term's hyperparameters,
  as a named list of single whole numbers. Empty, the default, leaves
  the number to the fitting layer. See
  [`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md).

- values:

  The values a path visits, for each hyperparameter the caller wrote
  out, as a named list of numeric vectors. Empty, the default, has the
  path build them. See
  [`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md).

- min_ratio:

  How far down the path over the size of the kink reaches, as a fraction
  of the value that empties the block: one number in \\(0, 1)\\, or
  `numeric(0)` for the fitting layer's own. See
  [`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md).

- search:

  How the term's own hyperparameters are covered when it has more than
  one carrying a kink: `"grid"` for every combination of them,
  `"cyclic"` for one at a time, or `character(0)` for the default. See
  [`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md).

- X:

  The design block, one row per observation: a numeric matrix or a
  two-dimensional Matrix. Empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- coef_names:

  The block's coefficient names, one per column, prefixed by `label`
  when there is one. Filled by
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md).

- blueprint:

  A named list of everything needed to reproduce the mapping on new
  rows. Its contents are the subclass's business; nothing outside the
  term reads them. Empty until
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  fills it.

- penalty:

  A penalties7 penalty on the block's coefficients, or `NULL` when the
  term is unpenalized.

- kind:

  One of `"seg"`, `"jump"` or `"jseg"`.

- var:

  The covariate expression, unevaluated.

- npsi:

  The number of break-points, an integer of at least 1.

- linear:

  Whether the block carries the linear effect.

- subformulas:

  A named list of one-sided formulas, one per developed coefficient.
  Empty where none is.

- spec:

  A named list of the resolved construction settings.

## Value

An S7 object of class `SegTerm`, inheriting from
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
with the six properties above beside the ten they supply.

## The six properties of its own

`kind` is `"seg"`, `"jump"` or `"jseg"`, and it decides almost
everything else: which coefficients the term has, whether its block is a
Jacobian, and which developments it will accept.

`var` is the covariate expression, `npsi` the number of break-points,
and `linear` whether the block carries the linear effect \\\beta x\\,
which `seg` and `jseg` do by default and `jump` never does.

`subformulas` holds one right-hand side per developed coefficient, and
`spec` the resolved construction settings: the starting positions, the
scaling factor `c0`, the restart budget `n_boot` and any smoother.

## Two kinds of block, and what distinguishes them

For
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
the block is the exact Jacobian: the contribution is differentiable in
\\\psi_k\\ away from the break-point, so the break-point is an ordinary
coefficient and a linear fit on the block is a Gauss-Newton step.
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
answers `TRUE`.

For
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
it is a working linearization with the weight \\W =
1/(2\lvert\tilde{x} - \psi\rvert)\\ frozen at the previous break-point,
and the break-point is **read off** two coefficients rather than
incremented.
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
answers `FALSE`, and a fitting layer routes such a term to exact working
fits alternating with committed read-offs. Smoothing the step
([`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)'s
`smoothed` argument) turns the break-point back into an ordinary
parameter, and the answer back to `TRUE`.

[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
is the one function that reports the position under either construction.

## See also

[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
for the three constructions;
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for the positions;
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for the block as they move;
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
for where to begin.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(120, 0, 10)))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)

# All three constructions are this class, and differ in `kind`.
vapply(list(seg(x), jump(x), jseg(x)), function(t) t@kind, character(1))
#> [1] "seg"  "jump" "jseg"

# Which decides the coefficients the term carries.
lapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
       function(t) term_coef_names(term_build(t, d)))
#> $seg
#> [1] "seg.beta"   "seg.gamma1" "seg.psi1"  
#> 
#> $jump
#> [1] "jump.delta1" "jump.g1"    
#> 
#> $jseg
#> [1] "jseg.beta"   "jseg.gamma1" "jseg.delta1" "jseg.g1"    
#> 

# And whether its block is a Jacobian or a working linearization.
vapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
       term_jacobian_block, logical(1))
#>   seg  jump  jseg 
#>  TRUE FALSE FALSE 
```
