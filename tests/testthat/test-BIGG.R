library(massdatabase)
x <- request_bigg_version()

test_that("request_bigg_version", {
  expect_equal(class(x), "data.frame")
  expect_equal(nrow(x), 3)
  expect_equal(ncol(x), 2)
})

y <- request_bigg_model_info()
test_that("request_bigg_model_info", {
  expect_equal(class(y), "data.frame")
  expect_equal(ncol(y), 5)
  expect_less_than(nrow(y), 200)
})

z <- request_bigg_universal_metabolite_info()
test_that("request_bigg_universal_metabolite_info", {
  expect_equal(class(z), "data.frame")
  expect_equal(ncol(z), 3)
  expect_less_than(nrow(z), 10000)
})


w1 <-
  request_bigg_universal_metabolite(metabolite_id = "g3p", return_form = "list")
w2 <-
  request_bigg_universal_metabolite(metabolite_id = "g3p", return_form = "data.frame")

test_that("request_bigg_universal_metabolite", {
  expect_equal(class(w1), "list")
  expect_equal(length(w1), 7)
  expect_equal(ncol(w2), 18)
  expect_error(request_bigg_universal_metabolite(metabolite_id = "g3p_test", return_form = "list"))
})
