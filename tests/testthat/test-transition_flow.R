test_that("compute_state_cost returns a valid cost matrix", {
  beta <- matrix(
    c(1, 0, 0, 0, 1, 0, 0, 0, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), paste0("g", 1:3))
  )

  cost <- compute_state_cost(beta)

  expect_equal(dim(cost), c(3, 3))
  expect_equal(unname(diag(cost)), c(0, 0, 0))
  expect_true(all(cost >= 0))
})

test_that("estimate_transition_flow matches requested marginals", {
  states <- c("resident", "inflammatory", "myofibroblast")
  source <- c(resident = 0.7, inflammatory = 0.2, myofibroblast = 0.1)
  target <- c(resident = 0.2, inflammatory = 0.4, myofibroblast = 0.4)
  cost <- matrix(
    c(0, 1, 2, 1, 0, 1, 2, 1, 0),
    nrow = 3,
    dimnames = list(states, states)
  )

  flow <- estimate_transition_flow(source, target, cost, lambda = 0.5)

  expect_true(flow$converged)
  expect_equal(rowSums(flow$flow), source, tolerance = 1e-6)
  expect_equal(colSums(flow$flow), target, tolerance = 1e-6)
  expect_true(is.finite(flow$expected_cost))
})

test_that("compute_fpi returns cell-level plasticity scores", {
  z <- matrix(
    c(0.8, 0.2, 0.1, 0.9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("c1", "c2"), c("resident", "inflammatory"))
  )
  flow <- matrix(
    c(0.5, 0.3, 0.1, 0.1),
    nrow = 2,
    dimnames = list(c("resident", "inflammatory"), c("resident", "inflammatory"))
  )

  fpi <- compute_fpi(z, flow = flow)

  expect_equal(nrow(fpi), 2)
  expect_true(all(c("entropy", "transition_potential", "fpi") %in% colnames(fpi)))
  expect_true(all(is.finite(fpi$fpi)))
})
