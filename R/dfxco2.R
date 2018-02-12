#' Detrended XCO2 seasonal cycles for GOSAT-ACOSv3.5 and DGVM simulations
#'
#' The detrended XCO2 seasonal cycle data presented herein, 
#' are the result of signal decomposition of XCO2 daily means by TransCom Region for 2009-2012.
#' Land fluxes (NBP) underwent simulated atmospheric transport,
#' using JAMSTEC's ACTM. XCO2 was then sampled via co-location to GOSAT observations.
#' Both GOSAT and simulated XCO2 underwent signal decomposition,
#' using the ccgcrv algoritm <https://www.esrl.noaa.gov/gmd/ccgg/mbl/crvfit/crvfit.html>.
#' 
#'
#' @docType data
#'
#' @usage data(dfxco2)
#'
#' @format A data.frame with 100136 rows and 16 variables:
#' \itemize{
#'   \item year: year (integer)
#'   \item month: month (integer)
#'   \item day: day (integer)
#'   \item func: values of the full function (harmonic+trend+short-term & long-term digital filter) (float)
#'   \item poly: values of the polynomial part of the function (float)
#'   \item smooth: values of the short-term smoothed curve; function + short-term filter of residualas (float)
#'   \item trend: values of the trend curve; polynomial plus long-term filter of residuals (float)
#'   \item detrend: values of the original data points minus the trend curve (float)
#'   \item smcycle: values of the smoothed, detrended annual cycle; smooth - trend (float)
#'   \item harm: values of the harmonic part of the function (float)
#'   \item smres: values of the short-term smoothed residuals from the function (float)
#'   \item trres: values of the long-term smoothed residual from the function (float)
#'   \item trres: values of the growth rate; the first derivative of the trend curve (float)
#' }
#' @keywords datasets
#'
#' @references Calle, L., B. Poulter, & P.K. Patra (2018)
#'
#' @source
#'
#' @examples
#' data(dfxco2)
"dfxco2"
