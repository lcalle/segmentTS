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
#' @format A data.frame object
#'
#' @keywords datasets
#'
#' @references Calle, L., B. Poulter, & P.K. Patra (2018)
#'
#' @source
#'
#' @examples
#' data(dfxco2)
"dfxco2"
