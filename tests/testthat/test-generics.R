

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
                        package = NULL)
  spec <- Mine(label = "mine", blueprint = list())
  expect_match(paste(capture.output(print(spec)), collapse = " "),
               "specification")

  built <- Mine(label = "mine", blueprint = list(order = 1:3))
  out <- paste(capture.output(print(built)), collapse = " ")
  expect_match(out, "MineTerm")
  expect_match(out, "built")
})


test_that("a structural class declaring nothing inherits the blueprint", {
  # The branch's contract is a declaration and not a convention: blueprint is
  # on the abstract structural_term, so a class written outside the package
  # carries one without knowing about it, and an unfilled one is an empty
  # list rather than a missing property.
  Bare <- S7::new_class("BareStructural", parent = structural_term,
                        package = NULL)
  b <- Bare(label = "bare")
  expect_true("blueprint" %in% S7::prop_names(b))
  expect_identical(b@blueprint, list())
  expect_false(term_is_built(b))
  expect_match(paste(capture.output(print(b)), collapse = " "), "specification")

  # and filling it is what being built means, on any such class
  filled <- Bare(label = "bare", blueprint = list(order = 1:3))
  expect_true(term_is_built(filled))
})

test_that("the three shipped structural classes inherit rather than restate", {
  # Removing the redeclaration must not move the property: S7 keeps a
  # redeclared property in the parent's slot, so the constructors' argument
  # order is the same either way, and every call in the package is by name.
  for (cl in list(GasTerm, RegimeTerm, MarginalBreakTerm)) {
    o <- cl(label = "x")
    expect_identical(o@blueprint, list())
    expect_false(term_is_built(o))
    # blueprint sits in the PARENT's slot rather than last, where each
    # subclass used to put it. The position is read from the parent rather
    # than written out, so adding a property to model_term moves both and
    # the assertion goes on saying the same thing.
    par <- which(names(structural_term@properties) == "blueprint")
    expect_length(par, 1L)
    expect_identical(S7::prop_names(o)[par], "blueprint")
    expect_true(par < length(S7::prop_names(o)))
  }
})
