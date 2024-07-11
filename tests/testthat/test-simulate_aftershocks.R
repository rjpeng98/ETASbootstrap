data(VCI_simulated_background_earthquakes)
set.seed(1)
out<-simulate_aftershocks(parameters_target = c(0.2424, 0.0068, 0.9771, 1.2200,
                                           0.0033, 2.4778, 0.1718),
                      background_catalog = VCI_simulated_background_earthquakes,
                      time_begin_background = "2000/01/01",
                      longitude_limit = c(-131, -126.25),
                      latitude_limit = c(48, 50),
                      time_limit = c("2000/01/01", '2018/04/27'),
                      magnitude_sample = VCI_magnitude_sample,
                      magnitude_threshold = 3.5,
                      time_zone="GMT")
test_that("data_frame", {
  expect_equal(is.data.frame(out), TRUE)
  expect_equal(tolower(names(out)),
               c("date", "time", "longitude", "latitude", "magnitude"))
  expect_equal(is.character(out$date), TRUE)
  expect_equal(is.character(out$time), TRUE)
  expect_equal(is.numeric(out$longitude), TRUE)
  expect_equal(is.numeric(out$latitude), TRUE)
  expect_equal(is.numeric(out$magnitude), TRUE)
  expect_equal(anyNA(out), FALSE)
})
