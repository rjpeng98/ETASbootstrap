data(VCI_earthquakes_plus)
data(VCI_simulated_background_earthquakes)
out_expect<- VCI_simulated_background_earthquakes
set.seed(1)
out<-simulate_background_earthquakes(VCI_earthquakes_plus)
test_that("test simulate_background_earthquakes", {
  expect_equal(out,out_expect )
})
