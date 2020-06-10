test_that("check coefficient length", {
  blblm_fit <- blblm(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris, m = 10, B = 50)
  co <- coef(blblm_fit)
  expect_equal(length(co), 3)
})
test_that("return object in blblm class with parallel", {
  blblm_fit <- blblm(Sepal.Length ~ Sepal.Width, data = iris, m = 10, B = 50, core = 2)
  expect_s3_class(blblm_fit, "blblm")
})
test_that("return object in blblm class without parallel", {
  blblm_fit <- blblm(Sepal.Length ~ Sepal.Width, data = iris, m = 10, B = 50)
  expect_s3_class(blblm_fit, "blblm")
})