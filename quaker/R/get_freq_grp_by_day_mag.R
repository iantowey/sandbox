#' Transforms geojson seismic object to continency table day bay magnitude
#'
#' @param seismic.geojson.obj GeoJson object
#'
#' @return A \code{\link[tibble]{tibble}} with frequency of earthquakes grouped by day and magnitude
#' @export
#'
#' @importFrom dplyr "%>%" transmute select group_by summarise
#' @importFrom tidyr spread
#'
#'
#' @seealso
#'  \code{\link{get_seismic_data}},
#'  \code{\link{fit_quaker}},
#'  \code{\link{flatten_to_table}},
#'  \code{\link{get_freq_grp_by_country_mag}},
#'  \code{\link{plot}}
#'
#' @examples
#' data <- get_freq_grp_by_day_mag(get_seismic_data(timeFrame = 'PAST_WEEK', minMagnitude = '1'))
#' data
#' data <- get_freq_grp_by_day_mag(get_seismic_data(timeFrame = 'PAST_MONTH', minMagnitude = '2.5'))
#' data
#'
get_freq_grp_by_day_mag <- function(seismic.geojson.obj){

  data <- seismic.geojson.obj
  #get frequency of earthquake magnitude grouped by day
  df.magnitudeByDay = data$features$properties %>%
    transmute(
      magnitude=cut(mag,breaks=c(0:ceiling(max(mag,na.rm = TRUE))), na.rm=TRUE),
      day=format(as.POSIXct(time/1000, origin="1970-01-01"),"%Y-%m-%d")) %>%
    group_by(day,magnitude) %>%
    summarise(count=n()) %>%
    spread(key=magnitude, value=count, fill=0)
  df.magnitudeByDay$total = apply(df.magnitudeByDay[,-c(1)],1,'sum')

  return(df.magnitudeByDay)
}

if(getRversion() >= "2.15.1") utils::globalVariables(c('n','n()','.','count','country','day','mag','magnitude','place'))
