library(testthat)
library(modelterms7)
# the smoother constructors a smooth term is built from. modelterms7 imports
# basis7 without attaching it, so a test writing bspline_smooth() bare needs
# it here; a user reaches it through library(statmodels7), which attaches
# every member.
library(basis7)

test_check("modelterms7")
