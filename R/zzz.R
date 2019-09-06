#' Example DataSet
#'
#' TMB.LBSPR EXAMPLE
#'
#' @format A ‘data.table’ with 25 obs. of  160 variables:
#' /describe{
#'   /item{Species}{Character string for species name}
#'   /item{Family}{Family-level taxonomic name. Important to use correct spelling if calling "StepwiseLH" for life history data.}
#'   /item{LH.Source}{Either "Study" or "Stepwise". Study=life history parameters inputted directly. Stepwise=Provide an Lmax value and obtained life history paramater esitmates through StepwiseLH package}
#'   /item{Par}{Name of the model parameter}
#'   /item{Val1}{First parameter of a probability distribution describing the corresponding parameter (e.g. mean)}
#'   /item{Val2}{Second parameter of a probability distribution describing the corresponding paramter (e.g. SD).}
#'   /item{Dist}{Name of the probability distribution to be used  (e.g. "Normal")}
#'   /item{Min}{Minimum bound for this parameter}
#'   /item{Max}{Maximum bound for this parameter}
#'   /item{Length_obs}{Lower edge of the size bins use to classify length observations}
#'   /item{Count_obs}{Counts or proportion of individual fish in a length bin. The following columns are simply bootstrap iterations of abundance-at-length data to incorporate uncertainty in length data into the model.}
#' }
#'
"Example"

#' @useDynLib TMB.LBSPR
.onUnload <- function (libpath) {
  library.dynam.unload("TMB.LBSPR", libpath)
}
