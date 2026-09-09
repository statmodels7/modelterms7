# The smoother constructors a smooth term is built from.
#
# modelterms7 Imports basis7 without attaching it, so a test writing
# bspline_smooth() bare cannot see it: an import is visible inside the
# package's namespace and not on the search path a test runs against. A user
# reaches these through library(statmodels7), which attaches every member.
#
# It is a helper rather than a line in tests/testthat.R because
# testthat::test_local() -- how the suite is run during the work -- sets up
# its own environment and never sources that file, so a line there is read
# only by R CMD check.
library(basis7)
