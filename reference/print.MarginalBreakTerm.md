# Print a Marginal Break-Point Term

Prints the label and the kind, how many latent break-points each group
carries, and, for a built term, over how many groups. A second line
lists the parameters, and a third names the prior where one was given.

## Arguments

- x:

  A
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md),
  built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. The method is called for the lines that it writes.

## Details

The form is

    <MarginalBreakTerm> 'jump' (jump): 1 latent break-point per group, integrated out (3 groups)
      parameters: psi1.mean, psi1.sd, delta1

The prior line appears only where `random(distrib = )` names a family;
under the default Gaussian there is nothing to name. A built term is
shown with its number of groups, and an unbuilt one is marked as a
specification.

## See also

[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)

# The first call prints a specification, and the second prints the same
# term built over three groups.
jump(x, psi ~ random(~ 1 | id), marginal = TRUE)
#> <MarginalBreakTerm> 'jump' (jump): 1 latent break-point per group, integrated out (specification)
#>   parameters: psi1.mean, psi1.sd, delta1
term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
#> <MarginalBreakTerm> 'jump' (jump): 1 latent break-point per group, integrated out (3 groups)
#>   parameters: psi1.mean, psi1.sd, delta1
```
