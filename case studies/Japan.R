# Install the ETAS and ETASbootstrap packages from CRAN
install.packages("ETAS")
install.packages("ETASbootstrap")

# Load the ETAS and ETASbootstrap packages
library(ETAS)
library(ETASbootstrap)

# Reshape the japan.quakes dataset
japan.quakes <- japan.quakes[, 1:5]
colnames(japan.quakes) <- c("date", "time", "longitude", "latitude", "magnitude")

ETAS_Boots(
  earthquake_data = japan.quakes,
  longitude_boundaries = c(128, 145),
  latitude_boundaries = c(27, 45),
  study_region = list(long = c(130, 135, 145, 140), lat = c(33, 30, 40, 43)),
  time_begin = "1926-01-08 00:00:00",
  study_start = "1953-05-26 00:00:00",
  study_end = "1990-01-08 00:00:00",
  magnitude_threshold = 5.5,
  time_zone = "GMT",
  round_off = TRUE,
  parameters_0 = NULL,
  number_simulations = 1000,
  confidence_level = 0.95,
  output_datasets = TRUE,
  output_estimates = TRUE
)