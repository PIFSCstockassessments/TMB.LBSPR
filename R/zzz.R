#' @useDynLib TMB.LBSPR
.onUnload <- function (libpath) {
  library.dynam.unload("TMB.LBSPR", libpath)
}
