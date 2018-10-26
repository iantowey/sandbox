#' Plot geojson object
#'
#' @param x GeoJson Object
#' @param ... other arguments to be handle internally by the plot function
#'
#' @return Plots an object of type seismic_geojson.
#' @export
#'
#' @importFrom dplyr "%>%" select
#' @importFrom tidyr spread
#' @importFrom leaflet leaflet addTiles setView clearShapes addCircles
#'
#' @seealso
#'  \code{\link{get_seismic_data}},
#'  \code{\link{fit_quaker}},
#'  \code{\link{flatten_to_table}},
#'  \code{\link{get_freq_grp_by_country_mag}},
#'  \code{\link{get_freq_grp_by_day_mag}}
#'
#' @examples
#' plot(get_seismic_data(timeFrame = 'PAST_WEEK', minMagnitude = '1'))
#' plot(get_seismic_data(timeFrame = 'PAST_MONTH', minMagnitude = '2.5'))
#'
plot.seismic_geojson = function(x, ...){

  x.flat <- flatten_to_table(x)

  leaflet(x.flat) %>%
    addTiles() %>%
    setView(lng = 0, lat = 0, zoom = 2) %>%
    clearShapes() %>%
    addCircles(
      lng=~longtitude,
      lat=~latitude,
      radius = ~10^magnitude/20,
      weight = 1,
      color = "#ff0000",
      fillOpacity = 0.5,
      popup = ~paste(
        "<ul>",
        "<li>Place:",place,"</li>",
        "<li>Magnitude:",magnitude,"</li>",
        "<li>Time:",as.POSIXct(time/1000, origin="1970-01-01"),"</li>",
        "<li>Depth:",depth," km</li>",
        "<li>Tsunami Risk:",ifelse(tsunami==1,"Y","N"),"</li>",
        "<li><a href=",url,">More Details...</></li>",
        "</ul>"
      ))
}

