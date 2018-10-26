#' Fit a Logisitc regression to the seismic data and anova analysis of fit
#'
#' @param seismic.geojson.obj GeoJson object
#'
#' @return Logisitc regression of earthquake data, check if longtitude, latitude, depth, magnitude are good predictors of tsunami occurance.
#' @export
#'
#' @importFrom stats anova binomial glm predict time
#' @importFrom dplyr "%>%" select
#' @importFrom tidyr spread
#'
#' @seealso
#'  \code{\link{get_seismic_data}},
#'  \code{\link{flatten_to_table}},
#'  \code{\link{get_freq_grp_by_country_mag}},
#'  \code{\link{get_freq_grp_by_day_mag}},
#'  \code{\link{plot}}
#' @examples
#'
#' data <- fit_quaker(get_seismic_data(timeFrame = 'PAST_WEEK', minMagnitude = '1'))
#' data
#' data <- fit_quaker(get_seismic_data(timeFrame = 'PAST_MONTH', minMagnitude = '2.5'))
#' data
fit_quaker <- function(seismic.geojson.obj){

  data <- seismic.geojson.obj
  dataset <- data$features %>% (function(x){
    mat <- do.call(rbind,data$features$geometry$coordinates)
    mat <- cbind(mat,data$features$properties$mag)
    mat <- cbind(mat,data$features$properties$tsunami)
    df <- as.data.frame(mat)
    colnames(df) = c('longtiude','latitude','depth','magnitude','tsunami')
    df
  })

  fit<-stats::glm(formula=tsunami~., data=dataset, family=stats::binomial(link='logit'))
  fit.anova <- stats::anova(fit, test="Chisq")

  list(
    fit=fit,
    fit.anova=fit.anova
    )
}
