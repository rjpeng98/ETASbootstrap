# Reshape the italy.quakes dataset
italy.quakes <- italy.quakes[, 1:5]
colnames(italy.quakes) <- c("date", "time", "longitude", "latitude", "magnitude")

# Italy I
ETAS_Boots(
  earthquake_data = italy.quakes,
  longitude_boundaries = c(6.15, 19),
  latitude_boundaries = c(35, 48),
  study_region = list(long = c(13, 16, 12, 9), lat = c(40, 42, 46, 44)),
  time_begin = "2005-04-16 00:00:00",
  study_start = "2009-03-01 00:00:00",
  study_end = "2013-11-01 00:00:00",
  magnitude_threshold = 3,
  time_zone = "GMT",
  round_off = TRUE,
  parameters_0 = NULL,
  number_simulations = 1000,
  confidence_level = 0.95,
  output_datasets = TRUE,
  output_estimates = TRUE
)

#Italy II
ETAS_Boots(
  earthquake_data = italy.quakes,
  longitude_boundaries = c(6.15, 19),
  latitude_boundaries = c(35, 48),
  study_region = list(long = c(9, 17.65, 17.65), lat = c(37, 37, 43)),
  time_begin = "2005-04-16 00:00:00",
  study_start = "2009-03-01 00:00:00",
  study_end = "2013-11-01 00:00:00",
  magnitude_threshold = 3,
  time_zone = "GMT",
  round_off = TRUE,
  parameters_0 = NULL,
  number_simulations = 1000,
  confidence_level = 0.95,
  output_datasets = TRUE,
  output_estimates = TRUE
)

