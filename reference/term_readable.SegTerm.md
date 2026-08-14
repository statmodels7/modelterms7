# What a Fitted Break-Point Term Is About

The quantities of the model the term defines, rather than the
coefficients of the working block it is fitted through: the linear
effect \\\beta\\ where the term carries one, the changes of slope
\\\gamma_k\\ and of level \\\delta_k\\, and the break-points \\\psi_k\\,
with the Jacobian from the coefficients so that a caller holding their
variance matrix can carry it across.

## Arguments

- term:

  A built
  [`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- zeta:

  The term's coefficients, in the order of its block.

- ...:

  Unused.

## Value

A list with `name`, `value`, `jacobian` and `scale`, as
[`term_readable`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
documents, or `NULL` where a coefficient carries a development.

## Details

A continuous term holds every one of those as a coefficient and the
Jacobian is a selection. A discontinuous one does not hold the
break-point at all: it is read off the auxiliary pair as \\\psi_k =
-g_k/\delta_k\\, which the joint construction reaches too, its quadratic
degenerating to that reading once the increment has vanished. The
Jacobian rows are then

\$\$\frac{\partial \psi_k}{\partial g_k} = -\frac{1}{\delta_k}, \qquad
\frac{\partial \psi_k}{\partial \delta_k} =
\frac{g_k}{\delta_k^{2}},\$\$

which is the delta method `segmented` reports a break-point's standard
error by. The joint construction reads its position from a quadratic
that also carries the change of slope, and the two readings differ by
the increment the iteration has left: the VALUE reported is the one the
term holds, so that this and
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
cannot give two numbers for one quantity, and the Jacobian is the fixed
point's, which the quadratic degenerates to once the increment has
vanished.

Every quantity is on the identity scale: a change is unbounded and a
break-point is a position on the covariate's own scale, held inside the
interval between the 5th and the 95th percentile. Where a coefficient
carries a development there is no single number to report – a
break-point then has one value per observation – and the method returns
nothing, leaving a caller to report the coefficients themselves.
