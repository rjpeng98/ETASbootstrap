# Define the broader region
broader_region <- list(lat = c(33, 42, 42, 33),
                       long = c(61, 61, 42, 42))
region.win <- spatstat.geom::owin(poly = list(x = broader_region$long, y = broader_region$lat))

# Extract earthquakes within the broader region
iran.quakes_North <- iran.quakes[
  spatstat.geom::inside.owin(iran.quakes$long, iran.quakes$lat, region.win), ]

# Reshape the extracted dataset
iran.quakes_North <- iran.quakes_North[, 1:5]
colnames(iran.quakes_North) <- c("date", "time", "longitude", "latitude", "magnitude")

ETAS_Boots(
  earthquake_data = iran.quakes_North,
  longitude_boundaries = c(42, 61),
  latitude_boundaries = c(33, 42),
  study_region = list(lat = c(38, 37.8, 37, 35, 34, 34.2, 36, 36.5, 38, 38, 37, 38, 40, 41),
                      long = c(43, 46, 47.5, 48, 51, 53, 57, 60, 59, 58, 55, 50, 48, 43)),
  time_begin = "1973-01-01 00:00:00",
  study_start = "1990-06-01 00:00:00",
  study_end = "2010-06-01 00:00:00",
  magnitude_threshold = 4.5,
  time_zone = "GMT",
  parameters_0 = NULL,
  number_simulations = 1000,
  confidence_level = 0.99,
  round_off = FALSE,
  output_datasets = TRUE,
  output_estimates = TRUE
)