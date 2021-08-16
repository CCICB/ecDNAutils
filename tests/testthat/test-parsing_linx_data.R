test_that("get_ecDNA_chains works", {
  expected_df <- structure(list(clusterId = 274L, chainId = 8L), row.names = c(NA, -1L), class = "data.frame")
  linx_dir = system.file(package="utilitybeltlinx", "CHP212")

  expect_equal(
    get_ecDNA_chains(linx_dir, "CHP212", verbose = FALSE),
    expected_df
  )
})
