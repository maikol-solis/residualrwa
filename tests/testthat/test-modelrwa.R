test_that("Check class", {
  set.seed(1)
  X1 <- rnorm(100)
  X2 <- rnorm(100)
  X3 <- rnorm(100)

  Y <- plogis(X2^3+10*X1*X2)>0.5

  data <- as.data.frame(cbind(Y,X1,X2,X3))

  ex1 <- modelrwa(
    response.name = "Y" ,
    control = NULL,
    fixed = NULL,
    variables = c("X1","X2", "X3"),
    data = data,
    family = binomial,
    include.interactions = TRUE
  )

  expect_identical(class(ex1), "modelrwa")
})
