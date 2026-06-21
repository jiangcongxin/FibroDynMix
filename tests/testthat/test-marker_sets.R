test_that("get_fibrodynmix_markers returns built-in marker priors", {
  human <- get_fibrodynmix_markers(species = "human", context = "scar")
  mouse <- get_fibrodynmix_markers(species = "mouse", context = "skin")

  expect_s3_class(human, "FibroDynMixMarkerSet")
  expect_true(all(c("resident", "inflammatory", "myofibroblast", "ECM-remodeling") %in% names(human)))
  expect_true("COL1A1" %in% human$`ECM-remodeling`)
  expect_true("Col1a1" %in% mouse$`ECM-remodeling`)
  expect_equal(attr(human, "species"), "human")
  expect_equal(attr(human, "context"), "scar")
})
