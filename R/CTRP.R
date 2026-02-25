#' @title CTRP drug response covariances data
#'
#' @description
#' The empirical per-drug dose-dose covariance from the 
#' Cancer Therapeutics Response Portal (CTRP v2, 2015) database 
#' \url{https://portals.broadinstitute.org/ctrp}. Based on a subset of 
#' n=374 drugs profiled at the same p=5 concentrations: 0.002, 0.016, 0.130, 
#' 1.000 and 8.300uM. Using replicate viability measurements (cell lines), we 
#' computed one p×p covariance matrix S(d) per drug, yielding n= 374 matrices.
#'
#' @examples
#' # Load the covariance data from the CTRP dataset
#' data("CTRP", package = "moewishart")
#' head(CTRP)
#'
"CTRP"
