#============================================#
#            F U N C T I O N S               #
#     ::required for segment analysis::      #
#============================================#
#' @import stats
NULL
#============================================================================================
#          Functions for segment-based statistics of event-type signals
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
segmentTS.1matchsignal <- function(obs.evnt, sim.evnt, limit4match = 0){
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
  for(i in 1:n){
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
  #ls.evnt = list(obs.evnt, sim.evnt)
  #names(ls.evnt) = c('obs.evnt', 'sim.evnt')
  return(list(obs.evnt=obs.evnt, sim.evnt=sim.evnt))
} 
