#'@title Compute bootstrap confidence intervals
#'
#'@importFrom stats runif
#'
#'@description A number (1000 by default) of earthquake data catalogs are
#'  simulated by bootstrap and recorded. A 2-D spatial and temporal ETAS model is
#'  fitted to each bootstrap-simulated
#'  earthquake data catalog, and the corresponding parameter estimates are
#'  recorded, which provides an empirical distribution for each estimate.
#'  For a given confidence level \eqn{1-\alpha} (0.95 by default), bootstrap
#'  confidence intervals are built from the empirical \eqn{\alpha/2} (0.025) and
#'  \eqn{1 -\alpha/2} (0.975) quantiles
#'  of the distributions of estimates for the parameters
#'  \eqn{(A,c,\alpha,p,D,q,\gamma)} of the ETAS model.
#'
#'@details Ogata (1998) proposed the 2-D spatial and temporal ETAS model, which is
#'  now widely used to decluster earthquake catalogs and, to a lesser extent,
#'  make short-term forecasts. In the 2-D spatial and temporal ETAS model, the
#'  behavior of the point process for which
#'  \eqn{\{(t_i,x_i,y_i,m_i),i=1,\dots,n\}} is a partial realization is
#'  characterized by the conditional intensity function
#'  \deqn{\lambda_{\beta,\mathbf{\theta}}(t,x,y,m
#'  \mid H_t) = s_{\beta}(m)\lambda_{\mathbf{\theta}}(t,x,y \mid H_t),}
#'  where \eqn{\beta} and \eqn{\mathbf{\theta} = (\nu,A,\alpha,c,p,q,D,\gamma)}
#'  are the model parameters and \eqn{s_\beta} is the probability density function
#'  (pdf) associated with the distribution of earthquake magnitudes. It is
#'  assumed that the distribution of the magnitude of earthquakes is independent
#'  of the joint distribution of the occurrence time of earthquakes and the 2-D
#'  spatial location of their epicenters. It can be expressed, for arbitrary
#'  \eqn{\beta \in (0, \infty)}, as \deqn{s_{\beta}(m) = \beta \exp \{
#'  -\beta(m-m_0)\},} where \eqn{m} and \eqn{m_0} represent the magnitude of the
#'  earthquake and the magnitude threshold, respectively.
#'  Moreover, \eqn{\lambda_{\mathbf{\theta}}(t,x,y \mid H_t)} represents the rate of
#'  observation of earthquakes in time and space, given the information on
#'  earthquakes prior to time \eqn{t}. This rate is expressed as the sum of two
#'  terms, namely
#'  \deqn{\lambda_{\mathbf{\theta}}(t,x,y \mid H_t) = \mu(x,y) + \sum_{i:t_i<t}k(m_i)g(t-t_i)f(x-x_i,y-y_i \mid m_i)}
#'   with
#'  \deqn{\mu(x,y) = \nu u(x,y),}where \eqn{\nu \in (0, \infty)}.
#'  The term \eqn{\mu(x,y)} is usually called "background seismicity rate" and represents the rate at which earthquakes independently occur around longitude \eqn{x} and latitude \eqn{y}.
#'  The \eqn{i}th term of the summation in \eqn{\lambda_{\theta}}, namely
#'  \deqn{k(m_i)g(t-t_i)f(x-x_i,y-y_i \mid m_i),}
#'  represents the effect of the \eqn{i}th earthquake before time \eqn{t} on the occurrence rate of earthquakes that would occur at time \eqn{t}, with an epicenter around \eqn{(x,y)}. Thus,
#'  \deqn{\sum_{i:t_i<t}k(m_i)g(t-t_i)f(x-x_i,y-y_i \mid m_i)} describes the total effect of all the earthquakes that occurred prior to time \eqn{t}, on the rate at which earthquakes would occur with an epicenter around \eqn{(x, y)} at time \eqn{t}.
#'  The expressions for \eqn{k}, \eqn{g}, and \eqn{f} are discussed individually as follows. First,
#'  \deqn{ k(m) = Ae^{\alpha(m-m_0)},\quad m \geq m_0 ,} can be interpreted as the expected number of earthquakes triggered by a previous earthquake with magnitude \eqn{m}, where \eqn{A \in (0, \infty)} and \eqn{\alpha \in (0, \infty)}. Second, for all \eqn{t \in (t_i, \infty)},
#'  \deqn{g(t-t_i) = \frac{p-1}{c} \, \left (1+\frac{t-t_i}{c} \right )^{-p}} is the pdf for the occurrence time of an earthquake triggered by the \eqn{i}th earthquake in the catalog, which occurred at time \eqn{t_i}, where \eqn{c \in (0, \infty)} and \eqn{p \in (1, \infty)}. Third,
#'  \deqn{f(x-x_i,y-y_i \mid m_i) = \frac{q-1}{\pi De^{\gamma(m_i-m_0)}} \, \left\{ 1+\frac{(x-x_i)^2+(y-y_i)^2}{De^{\gamma(m_i-m_0)}} \right\}^{-q}}
#'  is the pdf for the occurrence location (epicenter) of an earthquake triggered by the \eqn{i}th earthquake in the catalog, which occurred with magnitude \eqn{m_i} and an epicenter at \eqn{(x_i, y_i)}, where \eqn{D \in (0, \infty)}, \eqn{\gamma \in (0, \infty)}, and \eqn{q \in (1, \infty)}.
#'  For more details, see the articles of Zhuang et al. (2002, 2004). 
#'
#'@param earthquake_data
#' An object of class "data.frame" containing the following 5 columns:
#'\itemize{
#' \item date: Occurrence date of earthquakes in the format "yyyy-mm-dd"
#' \item time: Occurrence time of earthquakes in the format "hh:mm:ss"
#' \item longitude: Longitude of the epicenter of earthquakes in decimal degrees
#' \item latitude: Latitude of the epicenter of earthquakes in decimal degrees
#' \item magnitude: Magnitude of earthquakes (Any type of magnitude is accepted as far as it is used consistently and thoroughly.)
#'}
#'See VCI_earthquakes for an example with a rectangular study region; for a more general, polygonal study region, see JPN_earthquakes.
#'Note that West longitude and South latitude
#'  values should be negative, whereas East longitude and North latitude
#'  values are positive.
#'@param longitude_boundaries A numerical vector of length 2 (long_min, long_max)
#'  with the longitude boundaries of a rectangular space window,
#'  for which the earthquake catalog data are contained in \bold{earthquake_data}. If NULL (at the beginning of the execution of the program),
#'  long_min and long_max will be set (by the program) to the minimum and maximum values of the longitudes
#'  of earthquakes in \bold{earthquake_data}. 
#'  Together with \bold{latitude_boundaries}, \bold{longitude_boundaries} defines a region aimed to take edge effects into account in the analyses. 
#'  This region includes the study region and approximately 20\% more space around the study region.
#'@param latitude_boundaries A numerical vector of length 2 (lat_min, lat_max)
#'  with the latitude boundaries of a rectangular space window, for which
#'  the earthquake catalog data are contained in \bold{earthquake_data}. If NULL, lat_min and lat_max will be set to the minimum and maximum values of the latitudes of
#'  earthquakes in \bold{earthquake_data}.
#'@param study_region A list with two components (lat, long) of equal length specifying the coordinates
#'  of the vertices of a polygonal study region. The vertices must be written in anticlockwise order. If NULL, study_region will
#'  be filled with boundaries defining a rectangular space window 20\% narrower than the space window built from the longitude_boundaries and latitude_boundaries, while keeping the same center.
#'@param time_begin A character string, in the date-time format (yyyy-mm-dd
#'  hh:mm:ss), which identifies the start of the time span in
#'  \bold{earthquake_data}. If NULL, \bold{time_begin} will be set to the date-time of the
#'  first event in \bold{earthquake_data}.
#'
#'@param study_start A character string, in the date-time format, which
#'  specifies the start of the study period. If NULL, \bold{study_start} will be set to
#'  the date-time corresponding to that of \bold{time_begin} plus 20\% of the length of the
#'  time span in \bold{earthquake_data}.
#'@param study_end A character string, in the date-time format, which specifies the
#'  end of the study period. If NULL, it will be set to the date-time of the last event in
#'  \bold{earthquake_data}.
#'  Note: \bold{study_end} coincides with the end of the time span in
#'  \bold{earthquake_data}.
#'@param magnitude_threshold A decimal number which specifies the threshold to be used for the
#'  magnitudes of earthquakes. Only earthquakes with a magnitude greater than or
#'  equal to \bold{magnitude_threshold} will be considered, while the model
#'  is being fitted. If NULL, the minimum magnitude
#'  calculated from the events in \bold{earthquake_data} will be used for
#'  \bold{magnitude_threshold}.
#'
#'@param time_zone A character string specifying the time zone in
#'  which the occurrence times of earthquakes were recorded.
#'  The default "GMT"is the UTC (Universal Time Coordinates).
#'@param round_off A logical flag indicating whether or not to account for round-off error in coordinates of epicenters. 
#'@param parameters_0 A decimal vector of size 8 \eqn{(\nu, A, c, \alpha, p, D, q, \gamma)} to be used as an initial solution for the
#' iterative maximum likelihood estimation of the ETAS model parameters.
#' In particular, the values of parameters \eqn{\nu},
#'  \eqn{A}, \eqn{c}, \eqn{\alpha}, \eqn{D}, and \eqn{\gamma} are positive,
#'  and those of
#'  \eqn{p} and \eqn{q} are strictly greater than 1.
#'  If NULL, the values recommended by Ogata (1998) will be used.
#'@param number_simulations A positive integer which stands for the number of
#'  requested bootstrap simulations. The default value is 1000.
#'@param confidence_level A decimal number in (0, 1) which specifies the
#'  confidence level associated with the bootstrap confidence intervals that are built for
#'  the ETAS model parameters, and saved as outputs.
#'  It is set to 0.95 by default.
#'@param output_datasets A logical flag indicating whether or not
#' the bootstrap-simulated earthquake data catalogs must be written in CSV
#' files. The default setting is FALSE.
#'@param output_estimates A logical flag indicating whether or not
#'  the maximum likelihood estimates of parameters from each
#'  bootstrap-simulated earthquake data catalog must be written in a CSV file.
#'  The default setting is FALSE.
#'
#'
#'@return A list of three components:
#'\itemize{
#'   \item MLE: A numerical vector recording the maximum likelihood estimates of the ETAS model parameters
#'        \eqn{(A,c,\alpha,p,D,q,\gamma)}.
#'   \item ASE: A numerical vector recording the corresponding asymptotic standard errors.
#'   \item BootstrapCI: A matrix recording the corresponding bootstrap confidence intervals
#'   for
#'   the \bold{confidence_level} entered as input and all the other input arguments
#'   of the \link{ETAS_Boots} function starting with \bold{earthquake_data}.
#' }
#'@return When \bold{output_datasets}=TRUE, the simulated earthquake data
#'  catalogs are
#'  written in "Boot_N.csv", where "N" denotes the number of bootstrap
#'  simulation runs.
#'@return When \bold{output_estimates}=TRUE, the maximum likelihood
#'estimates of parameters
#'  from each simulated earthquake data catalog are written in "estimates.csv".
#'
#'@export
#'
#'@references
#'  Dutilleul, P., Genest, C., Peng, R., 2024. Bootstrapping for parameter uncertainty
#'  in the space-time epidemic-type aftershock sequence model. Geophysical Journal 
#'  International 236, 1601–1608.
#'  
#'  Ogata, Y. (1998). Space-time point-process models for earthquake
#'  occurrences. Annals of the Institute of Statistical Mathematics 50(2),
#'  379–402.
#'
#'  Zhuang, J., Y. Ogata, and D. Vere-Jones (2002). Stochastic declustering
#'  of space-time earthquake occurrences. Journal of the
#'  American Statistical Association 97(458), 369–380.
#'
#'  Zhuang, J., Y. Ogata, and D. Vere-Jones (2004). Analyzing earthquake
#'  clustering features by using stochastic reconstruction. Journal of
#'  Geophysical Research: Solid Earth 109(B05301).
#'
#'
#'@examples\donttest{
#'set.seed(23)
#'ETAS_Boots(earthquake_data = VCI_earthquakes,
#'           longitude_boundaries = c(-131, -126.25),
#'           latitude_boundaries = c(48, 50),
#'           study_region = list(long = c(-130.5, -130.5, -126.75, -126.75),
#'                               lat = c(49.75, 48.25, 48.25, 49.75)),
#'           time_begin = "2000/01/01 00:00:00",
#'           study_start = "2008/04/27 00:00:00",
#'           study_end = "2018/04/27 00:00:00",
#'           magnitude_threshold = 4,
#'           time_zone = "GMT",
#'           parameters_0 = c(0.65, 0.24, 0.0068, 0.97, 1.22, 0.0033, 2.48, 0.17),
#'           number_simulations = 4,
#'           confidence_level = 0.95,
#'           output_datasets = FALSE,
#'           output_estimates = FALSE)}
#'           
#'@examples\donttest{
#'ETAS_Boots(
#'  earthquake_data=JPN_earthquakes,
#'  longitude_boundaries = c(128, 145),
#'  latitude_boundaries = c(27, 45),
#'  study_region = list(long=c(130,135,145,140),
#'                     lat=c(33,30,40,43)),
#'  time_begin = "1926-01-08",
#'  study_start = "1953-05-26",
#'  study_end = "1990-01-08",
#'  magnitude_threshold = 5.5,
#'  time_zone = "GMT",
#'  round_off = FALSE,
#'  parameters_0 = c(0.524813924, 0.09, 0.045215442, 1.970176559, 
#'                   1.249620329, 0.002110203, 1.910492169,1.763149113 ),
#'  number_simulations = 2,
#'  confidence_level = 0.95,
#'  output_datasets = FALSE,
#'  output_estimates = FALSE
#')
#'
#'}           

ETAS_Boots<- function(earthquake_data,
                     longitude_boundaries=NULL,
                     latitude_boundaries=NULL,
                     study_region = NULL, 
                     time_begin=NULL,
                     study_start=NULL,
                     study_end=NULL,
                     magnitude_threshold=NULL,
                     time_zone="GMT",
                     round_off = FALSE,
                     parameters_0=NULL,
                     number_simulations=1000,
                     confidence_level=0.95,
                     output_datasets=FALSE,
                     output_estimates=FALSE)
{
  E_data<- as.data.frame(earthquake_data)
  long_range<- longitude_boundaries
  lat_range<- latitude_boundaries
  study_region<- study_region
  t_begin<- time_begin
  s_start<- study_start
  s_end<- study_end
  m_threshold<- magnitude_threshold
  tz<- time_zone
  param<- parameters_0
  n_sim<- number_simulations
  c_level<- confidence_level
  output_datasets<- output_datasets
  output_estimates<- output_estimates

  #check the format of earthquake_data

  dnames <- tolower(names(E_data))
  vnames <- c("date", "time", "longitude", "latitude", "magnitude")
  if (!identical(vnames, dnames))
    stop(paste("argument", sQuote("earthquake_data"),
               "must be a data frame with column names",
               toString(sQuote(vnames)),"(in this order)"))

  if (any(is.na(E_data[, vnames])))
    stop(paste(sQuote(vnames), "must not contain NA values"))

  if (!is.numeric(E_data[,3])
      || !is.numeric(E_data[,4])
      || !is.numeric(E_data[,5]))
    stop("longitude, latitude and magnitude columns must be numerical")


  # extract spatial coordinates and magnitude

  x <- E_data[,3]  # longitude of earthquake epicenters
  y <- E_data[,4]   # latitude of earthquake epicenters
  m <- E_data[,5]   # magnitude of earthquakes

  #check the argument related to spatial region
  
  if (is.null(long_range)){
    long_range<- c(min(x),max(x))
  } 
  
  if (is.null(lat_range)){
    lat_range<- c(min(y),max(y))
  }
  
  if (min(x) <long_range[1]
      ||max(x) >long_range[2]
      ||min(y) <lat_range[1]
      ||max(y) >lat_range[2]){
    stop("at least one earthquake in the given earthquake_data with epicenter outside
    the rectangular window constructed by longitude_rang and latitude_range ")
  }
  
  if (is.null(study_region)){
    long_study_min <- stats::median(long_range) - 0.5*(long_range[2]-long_range[1])*0.9
    long_study_max <- stats::median(long_range) + 0.5*(long_range[2]-long_range[1])*0.9
    lat_study_min <- stats::median(lat_range) - 0.5*(lat_range[2]-lat_range[1])*0.9
    lat_study_max <- stats::median(lat_range) + 0.5*(lat_range[2]-lat_range[1])*0.9
    study_region <- list(long = c(long_study_min, long_study_max, long_study_max, long_study_min),
                        lat = c(lat_study_min, lat_study_min, lat_study_max, lat_study_max))
  }
  
  if(setequal(names(study_region),c("lat","long")) ==FALSE){
    stop("study_region would be a list with two components named lat and long ")
  }
  
  if(length(study_region$long)!=length(study_region$lat)){
    stop("lat and long must have equal length")
  }
  
  if (!is.list(study_region)) {
    stop("study_region would be a list with components lat and long of equal length
         specifying the coordinates of the vertices of a polygonal study region. 
         The vertices must be listed in anticlockwise order.")
  }
  
  long_min<- min(study_region$long)
  long_max<- max(study_region$long)
  lat_min<- min(study_region$lat)
  lat_max<-max(study_region$lat)
  
  if(long_min < min(long_range) |
     long_max > max(long_range) |
     lat_min < min(lat_range) |
     lat_max > max(lat_range)) {
    stop("The study space window must be inside the boundaries constructed by longitude_boundaries and latitude_boundaries. ")
  }
  
  # check the format of arguments related to time
  # extract date and time of events
  if (is.character(t_begin)!=TRUE){
    stop(paste("time_begin should be a character string"))
  }
  if (is.character(s_start)!=TRUE){
    stop(paste("study_begin should be a character string"))
  }
  if (is.character(s_end)!=TRUE){
    stop(paste("study_end should be a character string"))
  }
  dt <- as.POSIXlt.character(paste(E_data$date, E_data$time), tz=tz)
  if (sum(duplicated(dt))> 0)
  {
    dtidx<- which(diff(dt) == 0)
    for (i in dtidx)
    {
      dt[i + 1] <- dt[i] + as.difftime(1, units="secs")
    }
    E_data$date<- substr(dt,1,10)
    E_data$time<- substr(dt,12,19)
    warning(paste("more than one event has occurred simultaneously!",
                  "\ncheck events", toString(dtidx),
                  "\nduplicated times have been altered by one second"))
  }

  if (is.unsorted(dt))
  {
    warning(paste("events were not chronologically sorted:",
                  "they have been sorted in ascending order"))
    E_data <- E_data[order(dt), ]
  }

  if (is.null(t_begin))
    t_begin<- paste(E_data$date[1],
                    E_data$time[1])

  if (is.null(s_end))
    s_end<- paste(E_data$date[nrow(E_data)],
                  E_data$time[nrow(E_data)])

  if (is.null(s_start)){
     s_start_str<-as.POSIXlt(ETAS::date2day(s_end,start=t_begin,tz=tz)*
                               0.2*24*60*60,origin =t_begin,tz=tz)
     s_start<-paste(substr(s_start_str,1,10),
                    substr(s_start_str,12,19))

  }
  #check format of magnitude_threshold, number_simulations, and confidence_level

  if (is.null(m_threshold)){
    m_threshold<-min(E_data[,5])}
  else if (!is.numeric(m_threshold)){
    stop("magnitude_threshold have to be numerical") }

  if (n_sim%%1!=0 || n_sim<=0){
    stop("number_simulations should be a positive integer")
  }

  if (c_level <= 0 || c_level >= 1){
    stop("confidence_level should be a decimal in (0, 1)")
  }


## ------------------------------------------------------------------------------------------------
  #rename the columns of E_data
  date<- E_data[,1]
  time<- E_data[,2]
  long<- E_data[,3]
  lat<- E_data[,4]
  mag<- E_data[,5]
  E_data<- data.frame(date,time,long,lat,mag)
  
  region.win <- spatstat.geom::owin(poly = list(x = study_region$long, 
                                                y = study_region$lat))
  
  S_total_cat<- ETAS::catalog(E_data,
                                  time.begin= t_begin,
                                  study.start= s_start,
                                  study.end= s_end,
                                  region.poly = study_region,
                                  mag.threshold=m_threshold,
                                  tz=tz,
                                  roundoff = round_off)


## ------------------------------------------------------------------------------------------------
  S.fit<- ETAS::etas(S_total_cat, param0 =param)

## ------------------------------------------------------------------------------------------------
#sample distribution of magnitudes
#(earthquakes in the study space-time window)
  t<- ETAS::date2day(paste(E_data$date,E_data$time), start= t_begin, tz= tz)
  t_min<- ETAS::date2day(s_start, start= t_begin, tz=tz)
  t_max<- ETAS::date2day(s_end, start= t_begin, tz=tz)
  E_data_new<- E_data[t>=t_min& t<=t_max,]
  
  mag_sample<- subset(E_data_new,
                      mag >= m_threshold
                      &spatstat.geom::inside.owin(x=long,y=lat, w=region.win))$mag


## ------------------------------------------------------------------------------------------------
  parameter<- c(as.numeric(S.fit$param[-1]))


## ------------------------------------------------------------------------------------------------
  estimatesSDS<-c()


## ------------------------------------------------------------------------------------------------
#G_0=background catalog(in the target space-time window)+auxiliary catalog
#pb: Prob. of being the background events
#bwd: bandwidth

  for(n in 1:n_sim)
    {

      quakes<-c()

      quakes<- data.frame(substr(S_total_cat$longlat.coord$dt,1,10),
                 substr(S_total_cat$longlat.coord$dt,12,19),
                 as.numeric(S_total_cat$longlat.coord[,1]),
                 as.numeric(S_total_cat$longlat.coord[,2]),
                 as.numeric(S_total_cat$revents[,4]+m_threshold),
                 as.numeric(S.fit$bwd),
                 as.numeric(S.fit$pb))

      colnames(quakes)<- c("date", "time", "longitude",
                           "latitude", "magnitude", "bandwidth",
                           "probability")

      bg_cata<- simulate_background_earthquakes(quakes)


  # create the clustering catalogue
  # after_cata: offspring in the target space-time window
      after_cata<-c()
      after_cata<-simulate_aftershocks(parameters_target= parameter,
                                   background_catalog= bg_cata,
                                   time_begin_background= t_begin,
                                   longitude_limit= long_range,
                                   latitude_limit= lat_range,
                                   time_limit= c(t_begin, s_end),
                                   magnitude_sample= mag_sample,
                                   magnitude_threshold= m_threshold,
                                   time_zone=tz
                             )

  #Merge the background and the clustering catalog together
      simu.cata<- c()
      simu.cata<- as.data.frame(rbind(as.matrix(after_cata),
                                      as.matrix(bg_cata)))
      simu.cata$long<- as.numeric(simu.cata$longitude)
      simu.cata$lat<- as.numeric(simu.cata$latitude)
      simu.cata$mag<- as.numeric(simu.cata$magnitude)

  #Sort it up

      simu.cata<-
        simu.cata[order(ETAS::date2day(paste(simu.cata$date,simu.cata$time),
                                                start= t_begin, tz=tz)),]


      simu.cata<-subset(simu.cata,is.na(date)==0&is.na(time)==0
                        &is.na(long)==0&is.na(lat)==0&is.na(mag)==0)


      S_total_cat_sim<- ETAS::catalog(simu.cata,
                                      time.begin= t_begin,
                                      study.start= s_start,
                                      study.end= s_end,
                                      region.poly = study_region,
                                      mag.threshold=m_threshold,
                                      roundoff = round_off)




      S.fit_sim<- ETAS::etas(S_total_cat_sim, as.vector(S.fit$param))

      estimatesSDS<- rbind(estimatesSDS,S.fit_sim$param)

    #output
      if(output_datasets){
        colnames(simu.cata)<- c("date","time","longitude","latitude","magnitude")
        utils::write.csv(simu.cata[,1:5],paste("Boot_",n,".csv"))}

      if(output_estimates)
        utils::write.csv(estimatesSDS[,-1],"estimates.csv")
  }
  MLE<- S.fit$param[-1]
  ASE<- S.fit$asd[S.fit$itr,-1]
  BootstrapCI <- cbind(stats::quantile(estimatesSDS[,2],
                                       c((1-c_level)/2,1-(1-c_level)/2)),
                   stats::quantile(estimatesSDS[,3],
                                       c((1-c_level)/2,1-(1-c_level)/2)),
                   stats::quantile(estimatesSDS[,4],
                                       c((1-c_level)/2,1-(1-c_level)/2)),
                   stats::quantile(estimatesSDS[,5],
                                       c((1-c_level)/2,1-(1-c_level)/2)),
                   stats::quantile(estimatesSDS[,6],
                                       c((1-c_level)/2,1-(1-c_level)/2)),
                   stats::quantile(estimatesSDS[,7],
                                       c((1-c_level)/2,1-(1-c_level)/2)),
                   stats::quantile(estimatesSDS[,8],
                                       c((1-c_level)/2,1-(1-c_level)/2)))
  colnames(BootstrapCI)<-c("A","c","alpha","p","D","q","gamma")
  out<-list(MLE,ASE,BootstrapCI)
  names(out)<-c("MLE","ASE","BootstrapCI")
  return(out)
}
