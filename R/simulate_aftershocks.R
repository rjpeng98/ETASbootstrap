#' @title Simulate a catalog of aftershocks
#'
#' @description When a catalog of background earthquakes is given, this function can be applied
#'   to simulate aftershocks under the intensity function
#'   \eqn{\sum_{i:t_i<t}\hat{k}(m_i)\hat{g}(t-t_i)\hat{f}(x-x_i,y-y_i \mid m_i)}, which
#'   is determined by the target parameter values given by the user.
#'
#' @param parameters_target A numerical vector of size 7, (\eqn{\hat{A},\hat{c},\hat{\alpha},\hat{p},\hat{D},\hat{q},\hat{\gamma}}), specifying the target values of
#'   parameters in the ETAS model.
#' @param background_catalog An object of class "data.frame" with 5 columns: recording date, time, longitude, latitude,
#'   and magnitude of the background events, in this order and in a format consistent with that of \bold{earthquake_data}
#'    in the function \code{\link{ETAS_Boots}}.
#' @param time_begin_background A character string, in the date-time format, that specifies the
#'   beginning of the time span in \bold{background_catalog}. If NULL, it will be set by the program to the
#'   date-time of the first earthquake in \bold{background_catalog}.
#' @param longitude_limit A vector of size 2 (xlim_min, xlim_max) specifying the
#'   longitude boundaries for the simulated aftershocks. If NULL, xlim_min and
#'   xlim_max will be set by the program to the minimum and maximum values of the
#'   longitude for the earthquakes in \bold{background_catalog}, respectively. Only the simulated
#'   aftershocks with a longitude inside \bold{longitude_limit} will be kept.
#' @param latitude_limit A vector of size 2 (ylim_min, ylim_max) specifying the
#'   latitude boundaries for the simulated aftershocks. If NULL, ylim_min and
#'   ylim_max will be set by the program to the minimum and maximum values of
#'   latitude for the earthquakes in \bold{background_catalog}, respectively. Only the simulated
#'   aftershocks with a latitude inside \bold{latitude_limit} will be kept.
#' @param time_limit A vector of size 2 (tlim_min, tlim_max) specifying the time span for
#'   the simulated aftershocks. If NULL, tlim_min and tlim_max will be set by the program to the
#'   date-time of the first and last earthquakes (in chronological order) in \bold{background_catalog},
#'   respectively. Only the simulated aftershocks inside the specified time span will be kept.
#' @param magnitude_sample A vector recording the sample from the distribution of earthquake magnitudes (\eqn{s_{\beta}(m)}).
#'   If NULL, the magnitudes of earthquakes in \bold{background_catalog} will be used.
#' @param magnitude_threshold A decimal value specifying the magnitude
#'   threshold to be applied. Only the simulated aftershocks with a magnitude of at least
#'   \bold{mag_threshold} will be kept. 
#'   If NULL, the minimum magnitude of the events in \bold{background_catalog} will be used as
#'   \bold{magnitude_threshold}.
#' @param time_zone A character string specifying the time zone. The default setting
#'   "GMT" is the UTC (Universal Time Coordinated).
#'
#' @return aftershocks_simulated: An object of class "data.frame" with 5 columns: recording the
#'   date, time, longitude, latitude and magnitude of the simulated
#'   aftershocks, in this order and a consistent format.
#' @export
#' 
#' @references 
#' Dutilleul, P., Genest, C., Peng, R., 2024. Bootstrapping for parameter uncertainty
#' in the space-time epidemic-type aftershock sequence model. Geophysical Journal 
#' International 236, 1601–1608. 
#'
#' @examples
#' set.seed(1)
#' simulate_aftershocks(parameters_target = c(0.2424, 0.0068, 0.9771, 1.2200, 
#'                                            0.0033, 2.4778, 0.1718),
#'                      background_catalog = VCI_simulated_background_earthquakes,
#'                      time_begin_background = "2000/01/01",
#'                      longitude_limit = c(-131, -126.25),
#'                      latitude_limit = c(48, 50),
#'                      time_limit = c("2000/01/01", "2018/04/27"),
#'                      magnitude_sample = VCI_magnitude_sample,
#'                      magnitude_threshold = 3.5,
#'                      time_zone="GMT")
#'


simulate_aftershocks<-function(parameters_target,
                               background_catalog,
                               time_begin_background=NULL,
                               longitude_limit=NULL,
                               latitude_limit=NULL,
                               time_limit=NULL,
                               magnitude_sample=NULL,
                               magnitude_threshold=NULL,
                               time_zone="GMT")
{
  param_target<- parameters_target
  bg_catalog<- background_catalog
  t_begin<- time_begin_background
  xlim<- longitude_limit
  ylim<- latitude_limit
  t_limit<- time_limit
  m_sample<- magnitude_sample
  m_threshold<- magnitude_threshold
  tz<- time_zone

  #check format of parameters_target
  if (!is.vector(param_target)
      || length(param_target)!=7
      || !is.numeric(param_target)){
    stop("parameters_target should be a numerical vector of size 7")
  }
  else if (param_target[1] <= 0
           || param_target[2] <= 0
           || param_target[3] <= 0
           || param_target[5] <= 0
           || param_target[7] <= 0
           || param_target[4] <= 1
           || param_target[6] <= 1){
    stop("values of parameters A, c, alpha, D, and gamma must be positive,
         and values of parameters p and q have to be greater than 1")
  }

  #check the format of background_catalog
  dnames <- tolower(names(bg_catalog))
  vnames <- c("date", "time", "longitude", "latitude", "magnitude")
  if (!identical(vnames,dnames))
    stop(paste("argument", sQuote("background_catalog"),
               "must be a data frame with column names",
               toString(sQuote(vnames)),"(in this order)"))

  if (any(is.na(bg_catalog[, vnames])))
    stop(paste(sQuote(vnames), "must not contain NA values"))

  if (!is.numeric(bg_catalog[,3])
      || !is.numeric(bg_catalog[,4])
      || !is.numeric(bg_catalog[,5]))
    stop("longitude, latitude and magnitude columns must be numerical")

  # extract date and time of events
  dt <- as.POSIXlt.character(paste(bg_catalog$date, bg_catalog$time), tz=tz)
  if (sum(duplicated(dt)) > 0)
  {
    dtidx<- which(diff(dt) == 0)
    for (i in dtidx)
    {
      dt[i + 1] <- dt[i] + as.difftime(1, units="secs")
    }
    bg_catalog$date<- substr(dt,1,10)
    bg_catalog$time<- substr(dt,12,19)
    warning(paste("more than one event has occurred simultaneously!",
                  "\ncheck events", toString(dtidx),
                  "\nduplicated times have been altered by one second"))
  }

  if (is.unsorted(dt))
  {
    warning(paste("events were not chronologically sorted:",
                  "they have been sorted in ascending order"))
    bg_catalog <- bg_catalog[order(dt), ]
  }

  #check format of time_begin_background
  if (is.null(t_begin)){
    t_begin<- paste(bg_catalog$date[1],bg_catalog$time[1])
  }else if (!is.character(t_begin)){
    stop("time_begin_background needs to be a character
         string format in the date-time format")
  }

  #check format of longitude_limit and latitude_limit
  if (is.null(xlim)){
    xlim<- c(min(bg_catalog$longitude), max(bg_catalog$longitude))
  }else if(!is.vector(xlim)
           || !is.numeric(xlim)
           || length(xlim)!=2
           || xlim[1]>=xlim[2]){
    stop("longitude_limit must be a vector of length 2 giving (xlim_min, xlim_max)")
  }

  if (is.null(ylim)){
    ylim<- c(min(bg_catalog$latitude), max(bg_catalog$latitude))
  }else if(!is.vector(ylim)
           || !is.numeric(ylim)
           || length(ylim)!=2
           || ylim[1]>=ylim[2]){
    stop("latitude_limit must be a vector of length 2 giving (ylim_min, ylim_max)")
  }

  #check format of time_limit
  if (is.null(t_limit)){
    t_limit<- c(paste(bg_catalog$date[1], bg_catalog$time[1]),
                paste(bg_catalog$date[nrow(bg_catalog)],
                      bg_catalog$time[nrow(bg_catalog)]))
  }else if(as.POSIXlt.character(t_limit[1]) >= as.POSIXlt.character(t_limit[2])
           ||!is.vector(t_limit)
           || length(t_limit)!=2){
    stop("time_limit must be a vector of size 2 (tlim_min,tlim_max),
         in which tlim_min and tlim_max should be character string
         in the data-time format")
  }

  #check format of magnitude_sample
  if (is.null(m_sample)){
    m_sample<- bg_catalog$magnitude
  }else if (!is.vector(m_sample)
            || !is.numeric(m_sample)){
    stop("magnitude_sample must be a numerical vector")
  }

  #check format of magnitude_threshold
  if (is.null(m_threshold)){
    m_threshold<-min(bg_catalog$magnitude)
    }else if (!is.numeric(m_threshold)){
    stop("magnitude_threshold have to be numerical")
    }

  #check format of time_zone
  if (is.null(tz)){
    tz<-"GMT"
  }else if (!is.character(tz)){
    stop("time zone have to be a character string")
  }


  para<-list(A=param_target[1],
             c=param_target[2],
             alpha=param_target[3],
             p=param_target[4],
             d=param_target[5],
             q=param_target[6],
             gamma=param_target[7])

  start.cata<- data.frame(as.numeric(bg_catalog[,3]),
                          as.numeric(bg_catalog[,4]),
                          as.numeric(bg_catalog[,5]),
                          ETAS::date2day(paste(bg_catalog[,1],
                                         bg_catalog[,2]),
                                         start=t_begin,
                                         tz=tz))

  #represent tlim in decimal days
  tlim<- c(ETAS::date2day(t_limit[1], start=t_begin, tz=tz),
           ETAS::date2day(t_limit[2], start=t_begin, tz=tz))

  rmag<- boot(m_sample)

  magthreshold<- m_threshold

  out<- simulate.etas(para, start.cata,
                 xlim, ylim, tlim, rmag, magthreshold)

  out<- as.data.frame(out[order(out[,4]),])

  names(out)[names(out) == "offspring.X"]<- "long"
  names(out)[names(out) == "offspring.Y"]<- "lat"
  names(out)[names(out) == "offspring.M"]<- "mag"
  names(out)[names(out) == "offspring.T"]<- "time"

  date_time<-as.POSIXlt(out$time*24*60*60,origin = t_begin, tz=tz)

  date<- substr(date_time,1,10)
  time<- substr(date_time,12,19)
  latitude<- out[,2]
  longitude<- out[,1]
  magnitude<- out[,3]

  after_simulate<- data.frame(date,time,longitude,latitude,magnitude)


  return(after_simulate)
}

#randomly allocating magnitudes to simulated aftershocks
boot<- function(x){
  function(n){ x[1 + stats::runif(n, 0, length(x))]}
}

#simulating aftershocks
simulate.etas<-function(para, start.cata,
                        xlim, ylim, tlim, rmag, magthreshold)
{

  initial.cata<- start.cata

  Noff<- stats::rpois(n=nrow(initial.cata),
                      lambda=para$A*exp(para$alpha*(start.cata[,3]-magthreshold)))
                      #Noff number of offsprings for each

  temp<-NULL

  Noff[is.null(Noff) | is.na(Noff)]<- 0


  #generate offsprings
  if(sum(Noff)>0){
    #print(paste(sum(Noff), "events generated."))

    dd<- rep(para$d*exp(para$gamma*(start.cata[,3]-magthreshold)),
             Noff)

    RR<- sqrt(dd*(stats::runif(sum(Noff))^(1/(1-para$q))-1))

    Theta<- stats::runif(sum(Noff))*pi*2

    offspring.X<- rep(start.cata[,1],Noff)
    offspring.X<- offspring.X  + RR*cos(Theta)/cos(mean(ylim)/180*pi)
    offspring.Y<- rep(start.cata[,2],Noff)
    offspring.Y<- offspring.Y  +RR*sin(Theta) #degree

    offspring.M<- rmag(sum(Noff))

    offspring.T<- rep(start.cata[,4],Noff)
    offspring.T<- offspring.T + (stats::runif(sum(Noff))^(1/(1-para$p))-1)*para$c


    temp<-data.frame(cbind(offspring.X, offspring.Y,
                           offspring.M, offspring.T))

    #select simulated earthquakes in the specified spcae-time window
    itag<- temp[,1]>=xlim[1] &temp[,1] <=xlim[2] &
           temp[,2]>=ylim[1] &temp[,2] <=ylim[2] &
           temp[,4]>=tlim[1] &temp[,4] <=tlim[2]

    if(sum(itag)>0){
      temp<- temp[itag,] #only keep events in the target space

      rm(Noff, offspring.X, offspring.Y, offspring.T, offspring.M) #remove
      temp1<- simulate.etas(para, temp,
                            xlim, ylim, tlim, rmag, magthreshold)
      temp<- data.frame (rbind(temp,temp1))

    }
    temp<- subset(temp,
                  temp[,1]>=xlim[1] &temp[,1] <=xlim[2] &
                  temp[,2]>=ylim[1] &temp[,2] <=ylim[2] &
                  temp[,4]>=tlim[1] &temp[,4] <=tlim[2])
  }

  temp
}

