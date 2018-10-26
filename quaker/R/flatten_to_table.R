#' flatten geojson object to data table with must salient covariates
#'
#' @param seismic.geojson.obj GeoJson object
#'
#' @return A flattened dataset to tabular format, useful for analysis.
#' @export
#'
#' @seealso
#'  \code{\link{get_seismic_data}},
#'  \code{\link{fit_quaker}},
#'  \code{\link{plot}},
#'  \code{\link{get_freq_grp_by_country_mag}},
#'  \code{\link{get_freq_grp_by_day_mag}}
#'
#' @examples
#' data <- flatten_to_table(get_seismic_data(timeFrame = 'PAST_DAY', minMagnitude = 'all'))
#' data
#' data <- flatten_to_table(get_seismic_data(timeFrame = 'PAST_WEEK', minMagnitude = '1'))
#' data
#' data <- flatten_to_table(get_seismic_data(timeFrame = 'PAST_MONTH', minMagnitude = '2.5'))
#' data
flatten_to_table <- function(seismic.geojson.obj){

  data <- seismic.geojson.obj

  #get summary quantiles of 'longtiude','latitude','depth','magnitude'
  resp =  data$features %>% (function(x){
    mat <- do.call(rbind,data$features$geometry$coordinates)
    df <- data.frame(mat)
    colnames(df) = c('longtitude','latitude','depth')
    df$magnitude <-  data$features$properties$mag
    df$tsunami <-  data$features$properties$tsunami
    df$time <-  data$features$properties$time
    df$place <-  data$features$properties$place
    df$url <-  data$features$properties$url
    df
  })

  return(resp)
}
