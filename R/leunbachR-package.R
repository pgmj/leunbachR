#' leunbachR: Leunbach Test Equating
#'
#' Implements the Leunbach test equating method, following the DIGRAM software
#' written by Svend Kreiner. Both direct and indirect equating are available,
#' with parametric bootstrap standard errors and diagnostic statistics
#' including the Goodman-Kruskal gamma test and orbit analysis for person fit.
#'
#' @keywords internal
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom graphics abline contour grid image lines par polygon
#' @importFrom stats complete.cases optimize pchisq pnorm quantile rmultinom setNames
#' @importFrom utils setTxtProgressBar txtProgressBar
"_PACKAGE"

# Silence R CMD check NOTEs for symbols injected into mirai::mirai() expressions.
# These names exist as parameters to the mirai task; the static analyzer can't
# see the binding.
utils::globalVariables(c("seed", "data_list"))
