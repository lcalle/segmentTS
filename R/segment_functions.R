#============================================#
#            F U N C T I O N S               #
#     ::required for segment analysis::      #
#============================================#
#' @import stats
NULL
#============================================================================================
# Functions for Series Distance of event-type signals (Ehret and Zehe 2011 HydroEarthSysSci)
# ..modified for full time series as in (Calle, Poulter, ana Patra 2018)
#============================================================================================

#' Match Signals
#'
#' This function matches the number of events in the time series. It is based
#' on Ehret and Zehe's (2011) conception of event-type signals. But for our purposes,
#' we treat the full time series as a single event and break up the time series into segements
#' categorized as separate signals.
#' This function is generally not used as we set all data points as the same 'event',
#' but we keep this function here for future application.
#' @param obs.evnt data.frame object with variables of start time (decimal.date), end time (decimal.date), match (integer)
#' @param sim.evnt data.frame, as above, but for simulated data
#' @param limit4match number to search the vector in obs.evnt,sim.evnt for similar events (set to zero)
#' @return a list with both obs.evnt,sim.evnt with an additional variable for the matching index; used in segmentTS.2catsignal
#' @export
segmentTS.1matchsignal     <- function(obs.evnt, sim.evnt, limit4match = 0){
  ##############################c
  #
  # D A T A   S T R U C T U R E
  # -- event data: i.e., peaks and troughs in the graph by visual inspection
  # -- (obs.evnt, sim.evnt) data.frames with variables:
  # --  ($ts) start time
  # --  ($te) end time
  # --  ($match) number of the matching observed event {obs$, sim$}
  #
  # Returns: for each 'o' in obs$ and 's' in sim$,
  #        : returns a list with both data.frames, incld. the number of matching event or -999 (no match)
  #
  ##############################c
  #determine overlap
  n=nrow(obs.evnt)
  m=nrow(sim.evnt)
  overlap <- matrix(data = -999, nrow = nrow(obs.evnt), ncol=nrow(sim.evnt))
  print(paste0("total obs = ",n))
  for(i in 1:n){
    if((i %% 50) == 0){Sys.sleep(0.001); print(paste0("i = ",i))}
    for(j in 1:m){
      overlap[i,j] = min(obs.evnt[i,'te'], sim.evnt[j,'te']) - max(obs.evnt[i,'ts'], sim.evnt[j,'ts']) + 1
      if(overlap[i,j] < limit4match){overlap[i,j]=-999}
    }
  }
  #end overlap
  
  #update objects (obs.evnt, sim.evnt) with the matching index
  while(max(overlap) > -999){
    indx = which(overlap==max(overlap), arr.ind=TRUE)
    obs.evnt$match[indx[1,1]] = indx[1,2]
    sim.evnt$match[indx[1,2]] = indx[1,1]
    overlap[indx[1,1],] = -999
    overlap[,indx[1,2]] = -999
  }
  ls.evnt = list(obs.evnt, sim.evnt)
  names(ls.evnt) = c('obs.evnt', 'sim.evnt')
  return(ls.evnt)
}

#' Categorize Signals
#'
#' This function takes in the time-series data for observed and simulated from segmentTS.1matchsignal.
#' Categorizes portions of the curve into clear signals {trough, up, no-event, down, or peak}
#' based on a difference equation to determine first and second derivatives.
#' @param dat data.frame object with variables of data value; variables derived from segmentTS.1matchsignal.
#' @param lolim lower limit for consideration of matching; set to low value so that all values potentially match.
#' @return data.frame object with variable 'pos', categorized as above; used in segmentTS.2eqsignal,
#' @export
segmentTS.2catsignal       <- function(dat, lolim = -999){  
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

#' Equalize Signals
#'
#' This function takes in the time-series data for observed and simulated from segmentTS.2catsignal.
#' Attempts to equalize the number of signals in simulated time-series to match number of signals in the obs.
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
#' @return list object with two data.frames, with number of peaks,troughs equalized; used in segmentTS.4segdist
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
  ls.eq = list(obs.eq, sim.eq)
  names(ls.eq) = c('obs.eq', 'sim.eq')
  return(ls.eq)
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
      ls.dat = list(o.peak, o.trough, s.peak, s.trough)
      return(ls.dat)
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
      
      #remove peak with lower value value
      if(!is.null(rm_val)){dfx <- dfx[!(dfx$val %in% rm_val),]}
      return(dfx)
}

#' Segment Distance
#'
#' This function takes in the time-series data for observed and simulated from segmentTS.3eqsignal.
#' Matches similar signals in the two time-series. Distance statistics are simulation - observation.
#' @param obs.seg data.frame object with variables derived from segmentTS.3eqsignal
#' @param sim.seg data.frame object, variables as in obs.seg, but for simulated data
#' @return list object with 4 outputs: time-series of the matching times and values (poly_t,poly) and the distance statistics (dist_tdiff,dist_vdiff)
#' @export
segmentTS.4segdist <- function(obs.seg, sim.seg){
  ##############################c
  #
  # D A T A   S T R U C T U R E
  # -- takes in the segmented data (i.e., points within each unique $pos factor on curve)
  # -- (obs.seg, sim.seg) segement time-series data
  # -- ($val) values, ordered by time
  # -- ($time) timing of values on curve, ordered by time
  #
  # Note: n= # of obs.segs; m= # of sim.segs
  # Returns: two dataframes; (1) dist.amp= (n) amplitude differences between obs and sims
  # -- focuses on main signals, not local maxima/minima
  ##############################c
  n=nrow(obs.seg)
  m=nrow(sim.seg)
  
  #create return vector
  dist_tdiff = rep(0,times = n)
  dist_vdiff = rep(0,times = n)
  
  #create equally spaced points within (sim.seg)
  poly_t.sim = seq.Date(from = sim.seg$time[1], to = sim.seg$time[m], length.out = n)
  poly_v.sim = approx(x= sim.seg$time, y= sim.seg$val, xout= poly_t.sim,  method='linear')$y

  poly_t.obs = seq.Date(from = obs.seg$time[1], to = obs.seg$time[n], length.out = n)
  poly_v.obs = approx(x= obs.seg$time, y= obs.seg$val, xout= poly_t.obs,  method='linear')$y

  #determine distance (error) btwn obs and sim curves/segments
  for(i in 1:n){
    dist_tdiff[i] = poly_t.sim[i] - poly_t.obs[i]
    dist_vdiff[i] = poly_v.sim[i] - poly_v.obs[i]
  }
  
  ls.diff = list(dist_tdiff, dist_vdiff, poly_t.sim, poly_v.sim)
  names(ls.diff) = c('dist_tdiff', 'dist_vdiff', 'poly_t', 'poly_v')
  return(ls.diff)
} 
