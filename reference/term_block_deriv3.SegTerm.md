# The Third and Fourth Derivatives of a Break-Point Term's Block

Return the third
([`term_block_deriv3()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv3.md))
and fourth
([`term_block_deriv4()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv4.md))
derivatives of the block of a smoothed
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
along three and four directions. A sharp construction returns zeros, for
the reasons that
[`term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.md)
records for its second derivative.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- coef:

  The coefficients, or `NULL` for the stored ones.

- v, u, s, r:

  Numeric vectors as long as the term's coefficients.

- ...:

  Unused.

## Value

A numeric matrix, one row per observation and one column per coefficient
of the term.

## Details

With \\u = x - \psi\\, the smoother \\s\\ of the term and \\d_j =
s^{(j)}/2\\, the columns of a smoothed term are \\H(u)\\ for a change
\\\gamma\\, \\S(u)\\ for a step \\\delta\\ and \\-(\gamma S(u) + \delta
d_2(u))\\ for a break-point, where \\H' = S = 1/2 + d_1\\. Write \\p_i\\
for the movement of the break-point along direction \\i\\ (its
development's design times the direction, gated where the position sits
against a confinement limit) and \\g_i\\, \\e_i\\ for the movements of
the change and of the step. A derivative in a break-point coefficient
multiplies a function of \\u\\ by \\-p_i\\ and raises its order by one,
and the change and the step enter the break-point column linearly. At
order three the columns are \$\$-d_3\\p_ap_bp_c,\qquad
-d_4\\p_ap_bp_c,\$\$ \$\$(\gamma d_4 + \delta d_5)\\p_ap_bp_c -
(g_ad_3 + e_ad_4)\\p_bp_c - (g_bd_3 + e_bd_4)\\p_ap_c - (g_cd_3 +
e_cd_4)\\p_ap_b,\$\$ and at order four \$\$d_4\\p_ap_bp_cp_e,\qquad
d_5\\p_ap_bp_cp_e,\$\$ \$\$-(\gamma d_5 + \delta d_6)\\p_ap_bp_cp_e +
\sum_i (g_id_4 + e_id_5) \prod\_{j\ne i}p_j.\$\$ The fourth derivative
of a term with a step reads the sixth derivative of the smoother, which
the smoothers of numericals7 do not carry; the method signals an error
there.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(40, 0, 10)))
b <- term_build(seg(x, smoothed = numericals7::smooth_probit()), d)
m <- length(term_coef_names(b))
dim(term_block_deriv3(b, v = rep(1, m), u = rep(1, m), s = rep(1, m)))
#> [1] 40  3
```
