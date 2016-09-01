#' Run synerdrug shiny interface
#'
#' @return Opens web browser on shiny interface.
#' @export
#'
#' @examples
#' \dontrun{
#' runSynerDrug()
#' }
runSynerDrug <- function(){
    shiny::runApp(system.file('SynerDrug',
                    package = 'synerdrug'))
}
