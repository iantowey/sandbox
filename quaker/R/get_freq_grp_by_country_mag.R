#' Transforms geojson seismic object to continency table day bay magnitude
#'
#' @param seismic.geojson.obj GeoJson object
#'
#' @return A \code{\link[tibble]{tibble}} with frequency of earthquakes grouped by day and
#' magnitude
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr "%>%" transmute select group_by summarise
#' @importFrom tidyr spread
#'
#'
#' @seealso
#'    \code{\link{get_seismic_data}},
#'    \code{\link{fit_quaker}},
#'    \code{\link{flatten_to_table}},
#'    \code{\link{get_freq_grp_by_day_mag}},
#'    \code{\link{plot}}
#'
#' @examples
#' data <- get_freq_grp_by_country_mag(
#'            get_seismic_data(timeFrame = 'PAST_WEEK', minMagnitude = '1')
#'         )
#' data
#' data <- get_freq_grp_by_country_mag(
#'            get_seismic_data(timeFrame = 'PAST_MONTH', minMagnitude = '2.5')
#'         )
#' data
#'
#' @export
get_freq_grp_by_country_mag <- function(seismic.geojson.obj){

  data <- seismic.geojson.obj

  if(is.null(states)){
    states = tolower(unname(unlist(fromJSON('https://gist.githubusercontent.com/mshafrir/2646763/raw/8b0dbb93521f5d6889502305335104218454c2bf/states_hash.json'))))
  }

  #get frequency of earthquake magnitude grouped by country, uses the char array 'states') array' to map earthquakes occuring in the USA to the country label 'USA'
  df.magnitudeByCountry = data$features$properties %>%
    transmute(
      magnitude=cut(mag,breaks=c(0:ceiling(max(mag,na.rm = TRUE))), na.rm=TRUE),
      country=place %>% strsplit(split=',') %>% sapply(., FUN = function(x){sub("^\\s+", "", x[2])}) %>% (function(x){ ifelse(tolower(x) %in% states,'USA',x) })) %>%
    group_by(country,magnitude) %>%
    summarise(count=n()) %>%
    spread(key=magnitude, value=count, fill = 0)
  df.magnitudeByCountry$total = apply(df.magnitudeByCountry[,-c(1)],1,'sum')

  return(df.magnitudeByCountry)
}
if(getRversion() >= "2.15.1") utils::globalVariables(c('n','n()','.','count','country','day','mag','magnitude','place'))
states = NULL
