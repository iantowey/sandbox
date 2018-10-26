#' Shiny app forthe quaker package
#'
#' @return Runs the shiny app
#' @export
#'
#' @examples
#' \dontrun{
#' quaker_shiny()
#' }
quaker_shiny = function() {
  appDir = system.file("shiny-examples", 'app', package = 'quaker')
  if (appDir == "") {
    stop("Could not find 'app' directory in package. Try re-installing `quaker`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}


