#============================================#
#            F U N C T I O N S               #
#     ::required for segment analysis::      #
#============================================#
#' @import stats
NULL
#============================================================================================
#          Functions for segment-based statistics of event-type signals
#============================================================================================

#' Equalize Signals
#'
#' This function takes in the time-series data for observed and simulated from segmentTS.2catsignal.
#' Attempts to equalize the number of peaks and troughs in simulated time-series,
#' to match number of signals in the observed time-series. This fn defines the boundary positions
#' for the segments in the full time-series.
#' Modify the criteria for removing false peaks/troughs in this section. 
#' i.e., remove short-term signals and focus on seasonal patterns of rise and fall.
#' Also, set individual criteria to smooth signal characterization for different simulated time-series.
#' @param obs.evnt data.frame object with variables derived from segmentTS.2catsignal
#' @param sim.evnt data.frame, variables as in obs.evnt, but for simulated data.
#' @param val.mindays integer number of timesteps (days) between peaks troughs; helps to remove false peaks and troughs.
#' @param manual_removal section of code identifying which peak or trough to remove from the time-series.
#' By visual inspection, if the first peak/trough is a false peak/trough,
#' then pass code as obs.peak[c(-1),] or obs.trough[c(-1),]. Add multiple removals via obs.peak[c(-1,-2,...),].
#' Pass code for multiple signals with semi-colons as in manual_removal="obs.peak[c(-1),]; sim.peak[c(-1),]"
#' @return list object with two data.frames with number of peaks,troughs equalized.
#' The data frame only contains the vector positions of the peaks and troughs.
#' @export
segmentTS.3eqsignal  <- function(obs.evnt, sim.evnt, val.mindays = 250, manual_removal = NULL){
  ###############################
  #
  # D A T A   S T R U C T U R E
  # -- (obs.evnt, sim.evnt) data.frame with variables:
  # -- ($pkval) peak values, ordered by time
  # -- ($tgval) trough values, ordered by time
  #
  # Note: n= # of obs.evnts; m= # of sim.evnts
  # Returns: list of two data.frames, with number of peaks,troughs equalized
  # -- focuses on main signals, not local maxima/minima
  ###############################
  #peak & trough values
  obs.peak   = obs.evnt[which(obs.evnt$pos== 2),]
  obs.trough = obs.evnt[which(obs.evnt$pos== -2),]
  sim.peak   = sim.evnt[which(sim.evnt$pos== 2),]
  sim.trough = sim.evnt[which(sim.evnt$pos== -2),]
 
  #-------------------------------------------
  # match segment by visual inspection
  # ..rise-rise, decline-decline
  # ..taken as input conditional code section
  #-------------------------------------------
  if(!is.null(manual_removal)){manual_removal}  

  #------------------------------------------
  # constrain peaks (+) troughs (-) values
  # & remove consecutive peaks, valleys
  #------------------------------------------
  #constrain peak values > 0 and trough values < 0
  #drop all rows where peak values < 0 || trough values > 0
  obs.peak   = obs.peak[which(obs.peak$val > 0),]
  obs.trough = obs.trough[which(obs.trough$val < 0),]
  sim.peak   = sim.peak[which(sim.peak$val > 0),]
  sim.trough = sim.trough[which(sim.trough$val < 0),]
  
  #remove false peaks and valleys
  #..keep only lowest of consective valleys and highest of consecutive peaks
  obs.eq = rm_false_peaksValleys(obs.peak,obs.trough)
  sim.eq = rm_false_peaksValleys(sim.peak,sim.trough)
  
  #peak & trough values
  obs.peak   = obs.eq[which(obs.eq$pos== 2),]
  obs.trough = obs.eq[which(obs.eq$pos== -2),]
  sim.peak   = sim.eq[which(sim.eq$pos== 2),]
  sim.trough = sim.eq[which(sim.eq$pos== -2),]
  
  #check for peak trough equality
  n=nrow(obs.peak)
  m=nrow(sim.peak)
  
  obs.seg= nrow(obs.peak) + nrow(obs.trough) - 1
  sim.seg= nrow(sim.peak) + nrow(sim.trough) - 1
  
  updateBool=FALSE
  if(obs.seg != sim.seg){
    ls.equal <- equalize_peaksValleys(o.peak   = obs.peak,
                                      o.trough = obs.trough,
                                      s.peak   = sim.peak,
                                      s.trough = sim.trough,
                                      type     = "minDays",
                                      min.days = val.mindays)
    updateBool=TRUE
  }
  
  #---------------------
  # format final data
  #---------------------
  if(updateBool==TRUE){
    #remove false peaks and valleys
    #..keep only lowest of consective valleys and highest of consecutive peaks
    obs.eq = rm_false_peaksValleys(ls.equal[[1]],ls.equal[[2]])
    sim.eq = rm_false_peaksValleys(ls.equal[[3]],ls.equal[[4]])
  }else{
    #reset data for removing false peaks
    obs.eq = rbind(obs.peak,obs.trough)
    sim.eq = rbind(sim.peak,sim.trough)
    #sort by date/time
    obs.eq = obs.eq[order(obs.eq$time),]
    sim.eq = sim.eq[order(sim.eq$time),]
  }
  
  #put data into list
  #ls.eq = list(obs.eq, sim.eq)
  #names(ls.eq) = c('obs.eq', 'sim.eq')
  return(list(obs.eq=obs.eq, sim.eq=sim.eq))
}
    # internal function used in SeriesDist.3EqualizeEvents
    rm_false_peaksValleys          <- function(df.peak,df.trough){
      
      #reset data for removing false peaks
      dfx = rbind(df.peak,df.trough)
      #sort by date/time
      dfx = dfx[order(dfx$time),]
      
      #keep only lowest of consective valleys and highest of consecutive peaks
      rm_vec = c()
      for(k in 1:(nrow(dfx)-1)){
        if(dfx$pos[k] == dfx$pos[k+1]){
          if(dfx$pos[k] == -2){
            #troughs: keep lowest value
            ifelse(dfx$val[k] < dfx$val[k+1], rm_vec <- c(rm_vec,-1*(k+1)), rm_vec <- c(rm_vec,-1*k))
          }else{
            #peaks: keep greatest value
            ifelse(dfx$val[k] > dfx$val[k+1], rm_vec <- c(rm_vec,-1*(k+1)), rm_vec <- c(rm_vec,-1*k))
          }
        }
      }
      #remove false peaks,valleys if they exist
      if(!is.null(rm_vec)){dfx = dfx[rm_vec,]}
      
      return(dfx)
}
    equalize_peaksValleys          <- function(o.peak, o.trough, s.peak, s.trough, type="equalPeaks or minDays",min.days=250){
      #o.peak is for observations
      #s.peak is for simulated data
      
      n= nrow(o.peak)
      m= nrow(s.peak)  
      #==========================
      # equalize peaks, troughs
      #==========================
      if(type == "equalPeaks"){
        for(i in 1:abs(n-m)){
          if((n-m) > 0){
            diff=1:(n-1)
            for(j in 1:(n-1)){
              diff[j]= (o.peak$val[j] - o.trough$val[j]) + (o.peak$val[j+1] - o.trough$val[j])
            } 
            #testing code
            j= which(min(diff) == diff,arr.ind = TRUE)
            #remove row with trough minimum value
            o.trough = o.trough[-j,]
            #remove row with peak minimum value
            rm_val = min(o.peak$val[j], o.peak$val[j+1])
            o.peak   = o.peak[which(o.peak$val != rm_val),]
            n=nrow(o.peak)
          }else if((n-m) < 0){
            diff=1:(m-1)
            for(j in 1:(m-1)){
              diff[j]= (s.peak$val[j] - s.trough$val[j]) + (s.peak$val[j+1] - s.trough$val[j])
            }
            j= which(min(diff) == diff,arr.ind = TRUE)
            #remove row with trough minimum value
            s.trough = s.trough[-j,]
            #remove row with peak minimum value
            rm_val = min(s.peak$val[j], s.peak$val[j+1])
            s.peak   = s.peak[which(s.peak$val != rm_val),]
            m=nrow(s.peak)
          }
        }
      }
      #====================================
      # min. period btwn periods, troughs
      #====================================
      if(type == "minDays"){
        o.peak   <- format_minDays_btwn_peakValley(o.peak,peakTrough = 'peak',min.days = min.days)
        o.trough <- format_minDays_btwn_peakValley(o.trough,peakTrough = 'trough',min.days = min.days)
        s.peak   <- format_minDays_btwn_peakValley(s.peak,peakTrough = 'peak',min.days = min.days)
        s.trough <- format_minDays_btwn_peakValley(s.trough,peakTrough = 'trough',min.days = min.days)
        
        n=nrow(o.peak)
        m=nrow(s.peak)
      }
      
      #store data for return
      #ls.dat = list(o.peak, o.trough, s.peak, s.trough)
      return(list(o.peak=o.peak, o.trough=o.trough, s.peak=s.peak, s.trough=s.trough))
}
    format_minDays_btwn_peakValley <- function(dfx, peakTrough = "peak or trough", min.days= 250){
      rm_val = NULL
      n= nrow(dfx)
      
      diff=1:(n-1)
      for(j in 1:(n-1)){
        diff[j]= abs(difftime(dfx$time[j],dfx$time[j+1], units='days'))
      } 
      
      j = which(diff < min.days,arr.ind = TRUE)
      if(length(j) > 0){
        rm_val=1:length(j)
        for(k in 1:length(j)){
          if(peakTrough == 'peak'){
            rm_val[k] = min(dfx$val[j[k]], dfx$val[j[k]+1])
          }else if(peakTrough == 'trough'){
            rm_val[k] = max(dfx$val[j[k]], dfx$val[j[k]+1])
          }
        }
      }#end do something if peaks/troughs closer than min.days
      
      #remove peak with lower value
      if(!is.null(rm_val)){dfx <- dfx[!(dfx$val %in% rm_val),]}
      return(dfx)
} 
