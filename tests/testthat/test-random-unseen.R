# A level the fit never saw: refused by default, a row of zeros on request,
# and the within-group design that a caller integrating its effect reads.

dd <- data.frame(x = c(0.5, 1, 2, -1, 3, 0.2),
                 g = factor(c("a", "b", "a", "c", "b", "c")))
b <- term_build(random(~ 1 + x | g), dd)
nd <- data.frame(x = c(3, 4, -2), g = c("new", "a", "other"))

test_that("an unseen level is refused unless asked for a zero row", {
  expect_error(term_predict(b, nd), "was not present at build time")
  Z <- as.matrix(term_predict(b, nd, unseen = "zero"))
  expect_identical(dim(Z), c(3L, term_npar(b)))
  expect_true(all(Z[c(1L, 3L), ] == 0))
  # a seen level gets the same row as without the option
  expect_identical(Z[2L, ], as.matrix(term_predict(b, nd[2L, ]))[1L, ])
})

test_that("the within-group design is the row one group's effect multiplies", {
  W <- term_within(b, nd)
  expect_identical(colnames(W), term_group(b)$names)
  expect_equal(unname(W), cbind(1, nd$x))
  # the block is the within design interacted with the indicators
  Z <- as.matrix(term_predict(b, dd))
  W0 <- term_within(b, dd)
  gi <- as.integer(dd$g)
  for (i in seq_len(nrow(dd))) {
    expect_equal(unname(Z[i, (gi[i] - 1L) * 2L + 1:2]), unname(W0[i, ]))
  }
  expect_null(term_within(term_build(linpar(~ x), dd), dd))
})
