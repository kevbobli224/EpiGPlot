#' Launch Shiny App For Package EpiGPlot
#'
#' A function that launches the shiny app for this package.
#' The shiny app allows us to easily import epigenetic factor data sets,
#' processing of data set if needed, and easily control plotting parameters
#' as well as plot features.
#'
#' @return No return value but open up a shiny page.
#'
#' @examples
#' \dontrun{
#' runEpiGPlot()
#' }
#'
#' @references
#' Medvedeva, Y. A., Lennartsson, A., Ehsani, R., Kulakovskiy, I. V., Vorontsov, I. E., Panahandeh, P., Khimulya, G., Kasukawa, T., &amp; Drabløs, F. (2015). Epifactors: A comprehensive database of human epigenetic factors and complexes. Database, 2015.
#' \href{http://database.oxfordjournals.org/content/2015/bav067.full}{Link}
#'
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp
runEpiGPlot <- function() {
    appDir <- system.file("shiny-scripts",
                          package = "EpiGPlot")
    shiny::runApp(appDir, display.mode = "normal")
    return()
}
