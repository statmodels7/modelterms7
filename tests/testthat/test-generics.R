

test_that("term_is_built answers for a structural term too", {
  d <- data.frame(y = rnorm(30), g = factor(rep(1:3, 10)))

  # An additive term is built when it has coefficient names; a structural one
  # contributes no design columns and so has none, and is built when its
  # blueprint is filled. The predicate asks each branch about its own.
  gs <- gas(p = 1, q = 1)
  gb <- term_build(gs, d)
  expect_false(term_is_built(gs))
  expect_true(term_is_built(gb))

  rs <- regime(k = 2)
  expect_false(term_is_built(rs))
  expect_true(term_is_built(term_build(rs, d)))

  # The additive branch is untouched, which is what makes this safe: every
  # .assert_built() caller is a method registered on an additive class.
  ls <- linpar(~ g)
  expect_false(term_is_built(ls))
  expect_true(term_is_built(term_build(ls, d)))

  expect_error(term_is_built("not a term"), "model_term")
})


test_that("the base print reports a structural class it does not know", {
  # Every shipped structural term registers a print method of its own, so
  # this line is reached only from outside the package -- and it must not
  # read `X`, which that branch has no property for.
  Mine <- S7::new_class("MineTerm", parent = structural_term,
                        properties = list(blueprint = S7::class_list),
                        package = NULL)
  spec <- Mine(label = "mine", blueprint = list())
  expect_match(paste(capture.output(print(spec)), collapse = " "),
               "specification")

  built <- Mine(label = "mine", blueprint = list(order = 1:3))
  out <- paste(capture.output(print(built)), collapse = " ")
  expect_match(out, "MineTerm")
  expect_match(out, "built")
})


test_that("a structural class carrying no blueprint answers rather than raises", {
  # `blueprint` is declared on additive_term and on each of the three shipped
  # structural classes, and NOT on the abstract structural_term, so the
  # branch's contract is a convention rather than a declaration. Reading the
  # property unconditionally raised S7's "Can't find property" from inside a
  # predicate whose own guard promises a logical for any model_term.
  Bare <- S7::new_class("BareStructural", parent = structural_term,
                        package = NULL)
  b <- Bare(label = "bare")
  expect_false(term_is_built(b))
  expect_match(paste(capture.output(print(b)), collapse = " "), "specification")
})
