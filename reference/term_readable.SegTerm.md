# The Quantities Reported for a Break-Point Term

Returns the quantities of the model that the term defines, in place of
the coefficients of the working block through which it is fitted: the
linear effect \\\beta\\ where the term carries one, the changes of slope
\\\gamma_k\\ and of level \\\delta_k\\, and the break-points \\\psi_k\\,
together with the Jacobian from the coefficients so that a caller
holding their variance matrix can carry it across.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- zeta:

  The term's coefficients, in the order of its block.

- ...:

  Unused.

## Value

A list with `name`, `value`, `jacobian` and `scale`, as
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
documents, or `NULL` if the list of quantities is empty.

## Details

A continuous or smoothed term holds every one of those as a coefficient,
and the Jacobian is a selection, as it is for a held term. A sharp
discontinuous term does not hold the break-point as a coefficient: it is
read off the auxiliary pair as \\\psi_k = -g_k/\delta_k\\, which the
joint construction reaches too, because its quadratic degenerates to
that reading once the increment has vanished. The Jacobian rows are then

\$\$\frac{\partial \psi_k}{\partial g_k} = -\frac{1}{\delta_k}, \qquad
\frac{\partial \psi_k}{\partial \delta_k} =
\frac{g_k}{\delta_k^{2}},\$\$

which is the delta method by which `segmented` reports a break-point's
standard error. The joint construction reads its position from a
quadratic that also carries the change of slope, and the two readings
differ by the increment that the iteration has left: the value reported
is the one that the term holds, so that this and
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
give the same number for one quantity, and the Jacobian is that of the
fixed point, to which the quadratic degenerates once the increment has
vanished.

Every quantity is on the identity scale: a change is unbounded and a
break-point is a position on the covariate's own scale, held inside the
interval between the 5th and the 95th percentile. Where a coefficient
carries a development there is no single number to report, because a
break-point then has one value per observation, and the method leaves
that coefficient to the caller to report. The exception is the position
of a held
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md):
its slot holds the development's coefficients \\p\\ in \\\psi_i =
w_i'p\\, and they are reported as `psi1.<column>`, one per column of the
sub-design.
