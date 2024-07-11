#' @title Simulate a catalog of background earthquakes
#'
#' @description In fitting the ETAS model to the earthquake data catalog of
#'   interest (\bold{earthquake_data}), the background intensity function \eqn{{\mu}}
#'   is estimated. This function performs a simulation
#'   of background events based on the estimate \eqn{\hat{{\mu}}}.
#'   The time period for the simulated background catalog is consistent with
#'   that of \bold{earthquake_data}.
#'
#' @param earthquake_data_plus An object of data.frame with 7 columns: date, time,
#'   longitude, latitude, magnitude, bandwidth, and probability, in this order and
#'   in a consistent format for the first 5 columns.
#'   The columns bandwidth and probability are two numeric vectors.
#'   The column bandwidth records the smoothness
#'   bandwidths used in variable kernel estimation and the column probability contains
#'   the probability for each earthquake in the catalog of interest (observed earthquakes)
#'   to be a background event; see the \link{etas} function in the ETAS package (Jalilian, 2019) and
#'   the articles of Zhuang et al. (2002, 2004).
#'
#' @return background_catalog: An object of data.frame with 5 columns: date,
#'   time, longitude, latitude, and magnitude of the simulated background earthquakes,
#'   in this order and a consistent format.
#'
#' @export
#'
#' @references 
#'  Dutilleul, P., Genest, C., Peng, R., 2024. Bootstrapping for parameter uncertainty
#'  in the space-time epidemic-type aftershock sequence model. Geophysical Journal 
#'  International 236, 1601–1608. 
#'  
#'  Jalilian, A. (2019). ETAS: An \R package for fitting the space-time ETAS model to earthquake data.
#'  Journal of Statistical Software, Code Snippets, 88(1), 1–39. doi:10.18637/jss.v088.c01.
#'  
#'  Zhuang, J., Y. Ogata, and D. Vere-Jones (2002). Stochastic declustering of
#'  space-time earthquake occurrences. Journal of the American Statistical
#'  Association 97(458), 369–380.
#'
#'  Zhuang, J., Y. Ogata, and D. Vere-Jones (2004). Analyzing earthquake
#'  clustering features by using stochastic reconstruction. Journal of
#'  Geophysical Research: Solid Earth 109(B05301).
#'
#' @examples
#' set.seed(1)
#' simulate_background_earthquakes(VCI_earthquakes_plus)
#'
#'
simulate_background_earthquakes<- function(earthquake_data_plus) {

  earthquake_data_plus<- as.data.frame(earthquake_data_plus)

  #check the format of argument VCI_earthquake_plus
  dnames <- tolower(names(earthquake_data_plus))

  vnames <- c("date",
              "time",
              "longitude",
              "latitude",
              "magnitude",
              "bandwidth",
              "probability")

  if (!identical(vnames,dnames))
    stop(paste("argument", sQuote("earthquake_data_plus"),
               "must be a data frame with column names",
               toString(sQuote(vnames)),"(in this order)"))


  if (any(is.na(earthquake_data_plus[, vnames]))){
    stop(paste(sQuote(vnames), "must not contain NA values"))}


  if (!is.numeric(earthquake_data_plus[,3]) ||
      !is.numeric(earthquake_data_plus[,4]) ||
      !is.numeric(earthquake_data_plus[,5]) ||
      !is.numeric(earthquake_data_plus[,6]) ||
      !is.numeric(earthquake_data_plus[,7]))
    stop("longitude, latitude, magnitude, bandwidth, and probability
         must be numerical")

  # extract date and time of events
  dt <- as.POSIXlt.character(paste(earthquake_data_plus$date,
                                   earthquake_data_plus$time))
  if (sum(duplicated(dt)) > 0)
  {
    dtidx<- which(diff(dt) == 0)
    for (i in dtidx)
    {
      dt[i + 1] <- dt[i] + as.difftime(1, units="secs")
    }
    earthquake_data_plus$date<- substr(dt,1,10)
    earthquake_data_plus$time<- substr(dt,12,19)
    warning(paste("more than one event has occurred simultaneously!",
                  "\ncheck events", toString(dtidx),
                  "\nduplicated times have been altered by one second"))
  }

  if (is.unsorted(dt))
  {
    warning(paste("events were not chronologically sorted:",
                  "they have been sorted in ascending order"))
    earthquake_data_plus <- earthquake_data_plus[order(dt), ]
  }


  #manipulating columns and plus one new column named Parents
  quakes<- data.frame(earthquake_data_plus[,3],
                      earthquake_data_plus[,4],
                      earthquake_data_plus[,5],
                      paste(earthquake_data_plus[,1],earthquake_data_plus[,2]),
                      earthquake_data_plus[,7],
                      1:nrow(earthquake_data_plus),
                      earthquake_data_plus[,6])

  colnames(quakes)<-c("long", "lat", "magnitude", "time", "Pb", "Parents", "bwd")

  #generate a background catalog
  #Parents==1 indicating the associated earthquake is an aftershock
  U<- c()
  U<- runif(length(quakes$long), 0, 1)
  quakes$Parents<- as.numeric(U >= quakes$Pb)

  #subset of background events
  G_0<- c()
  G_0<- subset(quakes, quakes[,6]==0)

  #add a Gaussian deviation to the location of background events
  G_d<- c()
  for(i in 1:length(G_0$bwd)){
    G_d<- rbind(G_d, MASS::mvrnorm(1,
                              mu= c(0,0),
                              Sigma= matrix(c(G_0$bwd[i]^2,0,
                                              0,G_0$bwd[i]^2),2)))
  }
  G_0$long<- G_0$long+G_d[,1]
  G_0$lat<- G_0$lat+G_d[,2]

  #reorder the locations of background events
  nrow<- sample(1:length(G_0$long), length(G_0$long), replace= F)
  G_0$long<-G_0$long[nrow]
  G_0$lat<-G_0$lat[nrow]

  #allocating new magnitude to simulated background events
  G_0$mag<- sample(quakes$magnitude, length(G_0$long), replace= T)

  G_0<-G_0[,1:4]

  #output
  out<- data.frame(substr(G_0$time,1,10),
                   substr(G_0$time,12,19),
                   as.numeric(G_0$long),
                   as.numeric(G_0$lat),
                   as.numeric(G_0$magnitude))
  colnames(out)<-c("date","time","longitude","latitude","magnitude")
  return(out)
}

