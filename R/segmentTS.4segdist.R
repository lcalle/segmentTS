#============================================#
#            F U N C T I O N S               #
#     ::required for segment analysis::      #
#============================================#
#' @import stats
NULL
#============================================================================================
#          Functions for segment-based statistics of event-type signals
#============================================================================================

#' Segment Distance
#'
#' This function takes in the full time-series data for observed and simulated from segmentTS.1matchsignal.
#' The boundaries of the segments in the time-series are defined using the
#' positions of the major peaks and troughs from fn segmentTS.3eqsignal(). Each segment
#' is then passed to segmentTS.4segdist to calculate the point-by-point segment statistics.
#' Matches similar signals in the two time-series. Distance statistics are simulation minus observation.
#' @param obs.seg data.frame object with variables derived from segmentTS.1matchsignal
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
  
  #ls.diff = list(dist_tdiff, dist_vdiff, poly_t.sim, poly_v.sim)
  #names(ls.diff) = c('dist_tdiff', 'dist_vdiff', 'poly_t', 'poly_v')
  return(list(dist_tdiff= dist_tdiff,
			  dist_vdiff= dist_vdiff,
			   poly_tsim= poly_tsim,
			   poly_vsim= poly_vsim,
			   poly_tobs= poly_tobs,
               poly_vobs= poly_vobs))
} 
