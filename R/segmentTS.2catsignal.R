#============================================#
#            F U N C T I O N S               #
#     ::required for segment analysis::      #
#============================================#
#' @import stats
NULL
#============================================================================================
#          Functions for segment-based statistics of event-type signals
#============================================================================================

#' Categorize Signals
#'
#' This function takes in the time-series data for observed and simulated from segmentTS.1matchsignal.
#' Categorizes portions of the curve into clear signals {trough, up, no-event, down, or peak}
#' based on a difference equation to determine first and second derivatives.
#' @param dat data.frame object with variables of data value; variables derived from segmentTS.1matchsignal.
#' @param lolim lower limit for consideration of matching; set to low value so that all values potentially match.
#' @return data.frame object with variable 'pos', categorized as above; used in segmentTS.2eqsignal,
#' @export
segmentTS.2catsignal <- function(dat, lolim = -999){  
  ##############################c
  #
  # D A T A   S T R U C T U R E
  # -- (in.data) time-series data (n = number of obs or simulated data points)
  # -- (in.data) data.frames with variables:
  # -- ($val) property value
  # -- ($pos) position on curve or 'case':
  #       trough{-2}; up{-1}; noEvent{0}; down{1}; peak{2}
  # -- lolim=minimum threshold for consideration in matching; set to -999
  #
  # Returns: data.frame with data for each 'k' in obs.data, returns the position of the curve in '$pos'
  #
  ##############################c
  
  #create 'pos' variable
  dat$pos = -3
  
  for(i in 2:(nrow(dat)-1)){
    if(dat$val[i] > lolim){
      #########main
      if( (dat$val[i] - dat$val[i-1]) < 0 & (dat$val[i+1] - dat$val[i]) > 0){
        dat$pos[i] = -2
      }else if( (dat$val[i] - dat$val[i-1]) < 0 & (dat$val[i+1] - dat$val[i]) < 0){
        dat$pos[i] = -1
      }else if( (dat$val[i] - dat$val[i-1]) > 0 & (dat$val[i+1] - dat$val[i]) > 0){
        dat$pos[i] = 1                
      }else if( (dat$val[i] - dat$val[i-1]) > 0 & (dat$val[i+1] - dat$val[i]) < 0){
        dat$pos[i] = 2
      }else{
        dat$pos[i] = 0     
      }
      #########end main
    }else{
      dat$pos[i] = 0
    }
  }#end forloop
  
  #set position on curve for first and last data point as the position determine on x=2 and x=n-1
  dat$pos[1]=dat$pos[2]
  dat$pos[nrow(dat)]=dat$pos[nrow(dat)-1]
  
  return(dat)
} 
