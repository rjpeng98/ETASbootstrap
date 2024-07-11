set.seed(23)
out1<- ETAS_Boots(earthquake_data = VCI_earthquakes,
           longitude_boundaries = c(-131, -126.25),
           latitude_boundaries = c(48, 50),
           study_region = list(long = c(-130.5, -130.5, -126.75, -126.75),
                               lat = c(49.75, 48.25, 48.25, 49.75)),
           time_begin = "2000/01/01",
           study_start = "2008/04/27",
           study_end = '2018/04/27',
           magnitude_threshold = 4,
           time_zone = "GMT",
           round_off = TRUE,
           parameters_0 = c(0.65, 0.24, 0.0068, 0.97, 1.22, 0.0033, 2.48, 0.17),
           number_simulations = 1,
           confidence_level = 0.95,
           output_datasets = FALSE,
           output_estimates = FALSE)
MLE<- c(0.085787129, 0.001166630, 1.608594123, 1.168492252,
        0.001965032, 1.967862690, 0.462813177 )
names(MLE)<- c("A","c","alpha","p","D","q","gamma")
ASE<- c(0.13559508, 0.25182557, 0.05161876, 0.01745401,
        0.32156358, 0.06711215, 0.27109897  )
names(ASE)<- c("A","c","alpha","p","D","q","gamma")
BootstrapCI <- cbind(c(0.02288696, 0.02288696),
                     c(0.001011814,0.001011814),
                     c(1.955256,1.955256),
                     c(1.173247,1.173247),
                     c(0.001562758,0.001562758),
                     c(2.956847,2.956847 ),
                     c(1.528977,1.528977))
colnames(BootstrapCI)<- c("A","c","alpha","p","D","q","gamma")
rownames(BootstrapCI)<- c("2.5%","97.5%")
out<-list(MLE,ASE,BootstrapCI)
names(out)<-c("MLE","ASE","BootstrapCI")
test_that("bootstrap confidence intervals", {
  expect_equal(out1$MLE,
               out$MLE,tolerance=0.1)
  expect_equal(out1$ASE,
               out$ASE,tolerance=0.1)
  expect_equal(out1$BootstrapCI,
               out$BootstrapCI,tolerance = 0.1)
})
