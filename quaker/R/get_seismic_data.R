#' Get geojson object from USGS
#'
#' @param timeFrame One of: PAST_HOUR, PAST_DAY, PAST_WEEK, PAST_MONTH
#' @param minMagnitude One of: 'all','1.0','2.5','4.5'
#'
#' @return A GeoJson object
#' @export
#'
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr "%>%" transmute select group_by summarize
#'
#' @seealso
#'  \code{\link{fit_quaker}},
#'  \code{\link{flatten_to_table}},
#'  \code{\link{get_freq_grp_by_country_mag}},
#'  \code{\link{get_freq_grp_by_day_mag}},
#'  \code{\link{plot}}
#' @examples
#'
#' data <- get_seismic_data(timeFrame = 'PAST_WEEK', minMagnitude = '1')
#' str(data)
#' data <- get_seismic_data(timeFrame = 'PAST_MONTH', minMagnitude = '2.5')
#' str(data)
get_seismic_data <- function(timeFrame = c('PAST_HOUR','PAST_DAY','PAST_WEEK','PAST_MONTH'),
                           minMagnitude = c('all','1.0','2.5','4.5')){
  #get time frame argument value
  timeFrame.arg = match.arg(timeFrame)

  #get minimum magnitude argument value
  minMagnitude.arg = match.arg(minMagnitude)

  #
  timeFrame.arg.normalized = timeFrame.arg %>% gsub(pattern='PAST_',replacement='',x=.) %>% tolower

  #build url
  url = paste0('http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/',minMagnitude.arg,'_',timeFrame.arg.normalized,'.geojson')

  #fetch geojsin from url
  seismic.geojson = tryCatch({
    fromJSON(url)},warning = function(w) {
      print(paste("WARNING",w))
    }, error = function(e) {
      print(paste("ERROR",e, url))
    }, finally = {
      NULL
    })

  class(seismic.geojson) <- 'seismic_geojson'
  return(seismic.geojson)
}
if(getRversion() >= "2.15.1") utils::globalVariables(c('n','n()','.','count','country','day','mag','magnitude','place'))
