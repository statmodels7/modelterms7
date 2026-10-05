# gas(scaling = d) stores the exponent of the expected information the score
# is multiplied by; the term itself only carries it, the fitting layer reads it.

test_that("gas() stores the scaling and validates it", {
  expect_identical(gas()@scaling, 0)
  expect_identical(gas(scaling = 1)@scaling, 1)
  expect_identical(gas(scaling = -0.25)@scaling, -0.25)
  expect_error(gas(scaling = "inverse"), "single finite number")
  expect_error(gas(scaling = Inf), "single finite number")
  expect_error(gas(scaling = numeric(0)), "single finite number")
})

test_that("a scaled term says so when printed, an unscaled one does not", {
  expect_output(print(gas(scaling = 0.5)), "information to the power -0.5")
  out <- utils::capture.output(print(gas()))
  expect_false(any(grepl("information", out)))
})
