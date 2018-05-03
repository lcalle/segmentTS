#============================================#
#            F U N C T I O N S               #
#     ::required for segment analysis::      #
#============================================#
#' @import stats
NULL

#' Segment-based Statistics & Figures
#'
#' This function collates segment-based statistics based on the 
#' observed and simulated time-series.
#' @import graphics
#' @import grDevices
#' @param obs.evnt data.frame object for full time-series of observational data, from segmentTS.1matchsignal.
#' @param sim.evnt data.frame object for full time-series of simulated data, from segmentTS.1matchsignal.
#' @param ls.evnt.pos list object from segmentTS.3eqsignal with positions of
#' the peaks and troughs in the full time-series. Defines the boundaries of 
#' the segments.
#' @param obs.name name of observational data for table (string).
#' @param sim.name name of simulated data for table (string).
#' @param time.units units of time for calculating period length in the time-series; default is days.
#' @param val.units units of the value for the variable in the time-series.
#' @param save.plot save plot as a pdf (TRUE/FALSE); default is FALSE; default out is the current working directory.
#' @param outDir location to save plot; default location is getwd()).
#' @param region.name name of underlying region for table (string).
#' Default is 'null_region'; not important unless evaluating multiple regions.
#' At minimum, requires variables of values and date (YYYY-MM-DD)
#' @return a data.frame object segment-based statistical summaries for all segments, obs and sim.
#' @export
segmentTS.4statsplots <- function(obs.evnt,sim.evnt,ls.evnt.pos,obs.name="obs",sim.name="sim",time.units='days',val.units=NULL,save.plot=FALSE,outDir=getwd(),region.name=NULL){

	#separate observation,simulation
	obs.seg.time = ls.evnt.pos[['obs.eq']]
	sim.seg.time = ls.evnt.pos[['sim.eq']]

	#----------------------------------------
	# Get number of segments and Store Data
	#----------------------------------------
	num_segments_obs = nrow(obs.seg.time)-1
	num_segments_sim = nrow(sim.seg.time)-1

	num_segments = min(num_segments_obs, num_segments_sim)
    
    #empty data frame 
    df_errorStats= data.frame( region  = rep(ifelse(!is.null(region.name),region.name,'null_region'), times=num_segments*2),
                               tracer  = rep(c(obs.name,sim.name),times=num_segments),
                               n_obs   = rep(0,num_segments*2),
                               segment = rep(1:num_segments, each=2),
							   segment_type = rep(NA,num_segments*2),
							   period           = rep(0,num_segments*2),
                               amplitude        = rep(0,num_segments*2),
                               timing_rmse         = rep(0,num_segments*2),
                               timing_bias_max     = rep(0,num_segments*2),
                               timing_bias_min     = rep(0,num_segments*2),
                               timing_bias_mean    = rep(0,num_segments*2),
                               timing_bias_var     = rep(0,num_segments*2),
                               magnitude_rmse      = rep(0,num_segments*2),
                               magnitude_bias_max  = rep(0,num_segments*2),
                               magnitude_bias_min  = rep(0,num_segments*2),
                               magnitude_bias_mean = rep(0,num_segments*2),
                               magnitude_bias_var  = rep(0,num_segments*2))
   for(seg_j in 1:num_segments){
      seg_beg = seg_j
      seg_end = seg_j+1
      
      #-------------------------------------------------------
      # get segment based on boundary positions for segment
      #-------------------------------------------------------
      obs.segment = obs.evnt[which(obs.evnt$time >= obs.seg.time$time[seg_beg] & obs.evnt$time <= obs.seg.time$time[seg_end]),]
      sim.segment = sim.evnt[which(sim.evnt$time >= sim.seg.time$time[seg_beg] & sim.evnt$time <= sim.seg.time$time[seg_end]),]
      
      #---------------------------------------
      # get segment-based distance statistics
      #---------------------------------------
      segment.dist = segmentTS.4segdist(obs.seg = obs.segment, sim.seg = sim.segment)
      
      #############################
      #          P L O T          #
      #  multi-plot evaluations   #
      #############################
      #tri-plot with errors for segment
      #..only create one plot per model for the seasonal cycle
      if(save.plot==TRUE){grDevices::pdf(paste0(outDir,"/plot_",sim.name,"_segmentTS.pdf"),height=10,width=18)}
      graphics::par(mar=c(4,6,5,3),mfcol=c(2,2))

      #plot attributes
      x.lim = c(min(obs.evnt$time[1],sim.evnt$time[1]), max(obs.evnt$time[length(obs.evnt$time)],sim.evnt$time[length(sim.evnt$time)]))
      y.lim = c(min(obs.evnt$val,sim.evnt$val), max(obs.evnt$val,sim.evnt$val))
      plot.dates = seq.Date(from = x.lim[1], to = x.lim[2], by = 'month')
      plot.all.dates = seq.Date(from = obs.evnt$time[1], to = obs.evnt$time[length(obs.evnt$time)], by = 'month')

      #plot
      graphics::plot(x= obs.evnt$time, y=  obs.evnt$val,lwd=3,
             type='l', col='black', ylab='XCO2 (ppm)', xlab='',
             xlim= x.lim,
             ylim= y.lim, xaxt='n', cex.lab=2,cex.axis=1.5,
             cex.main=2,
             main=paste0('observation (black), simulation (red)',ifelse(!is.null(region.name),paste0('\n',region.name,''))))
      graphics::axis(side= 1, at= plot.all.dates, labels= format(plot.all.dates, "%b %y"), cex.axis = 1.5)
      graphics::lines(x= sim.evnt$time,y= sim.evnt$val, type='l', col='red', lwd=2)

      #points for matching segments
      graphics::points(x=obs.segment$time[seq(1,length(obs.segment$time), length=25)], y=obs.segment$val[seq(1,length(obs.segment$time), length=25)], pch=21,cex=1.2, bg='deepskyblue')
      graphics::points(x=sim.segment$time[seq(1,length(sim.segment$time), length=25)], y=sim.segment$val[seq(1,length(sim.segment$time), length=25)], pch=21,cex=1.2, bg='deepskyblue')
      graphics::abline(h=0,lty=1, lwd=0.5, col='black')

      #----------------------------------------
      # plot Close-up of Matching Segments
      #----------------------------------------
      #plot info for clarity
      x.main = sim.evnt[which((sim.evnt$time %in% obs.segment$time)==TRUE),'time']
      y.main = sim.evnt[which((sim.evnt$time %in% obs.segment$time)==TRUE), 'val']
      x.lim = c(min(obs.segment$time[1],sim.segment$time[1]), max(obs.segment$time[length(obs.segment$time)],sim.segment$time[length(sim.segment$time)]))
      y.lim = c(min(obs.segment$val,sim.segment$val), max(obs.segment$val,sim.segment$val))
      main.label = 'Matching Segment'
      col.sim    = 'red'

      #close up plots
      graphics::plot(x= x.main, y= y.main,
             type='l', col=col.sim, ylab=ifelse(!is.null(val.units),val.units,'units'), xlab='',xlim= x.lim, ylim= y.lim, xaxt='n', cex.lab=2,cex.axis=1.5,
             cex.main=2,main=main.label)
      graphics::axis(side= 1, at= plot.dates, labels= format(plot.dates, "%b %y"), cex.axis = 1.5)
      graphics::lines(x= obs.evnt$time, y=obs.evnt$val,type='l')
      graphics::lines(x= obs.segment$time, y=obs.segment$val,type='l', col='black', lwd=4)
      graphics::lines(x= sim.segment$time, y=sim.segment$val,type='l', col=col.sim, lwd=4)
      graphics::points(x=obs.segment$time[seq(1,length(obs.segment$time), length=35)], y=obs.segment$val[seq(1,length(obs.segment$time), length=35)], pch=21,cex=1.2, bg='deepskyblue')
      graphics::points(x=sim.segment$time[seq(1,length(sim.segment$time), length=35)], y=sim.segment$val[seq(1,length(sim.segment$time), length=35)], pch=21,cex=1.2, bg='deepskyblue')
      #plot only a few connecting segments btwn obs and sim, otherwise clutter
      graphics::segments(x0= segment.dist[['poly_t.obs']][c(TRUE,rep(FALSE,5))], y0= segment.dist[['poly_v.obs']][c(TRUE,rep(FALSE,5))],
                 x1= segment.dist[['poly_t.sim']][c(TRUE,rep(FALSE,5))],
                 y1= segment.dist[['poly_v.sim']][c(TRUE,rep(FALSE,5))], col='grey75')
      graphics::abline(h=0,lty=1, lwd=0.5, col='black')

      #----------------------------------------
      # Plot Timing (phase) & Magnitude Error
      #----------------------------------------
      #plot phase errror
      graphics::plot(x=obs.segment$time, xaxt='n', y=segment.dist[['dist_tdiff']], type='o', pch=21, bg='grey25', ylab='Days', xlab='', xlim= x.lim,
             cex.lab=2, cex.axis=1.5, cex.main=2, main='Phase Eror')
      graphics::axis(side= 1, at= plot.dates, labels= format(plot.dates, "%b %y"), cex.axis = 1.5)
      graphics::abline(h=0,lty=1, lwd=0.5, col='black')
      #plot magnitude error
      graphics::plot(x=obs.segment$time, xaxt='n',  y=segment.dist[['dist_vdiff']], type='o', pch=24, bg='grey25', ylab='ppm', xlab='', xlim= x.lim,
             cex.lab=2, cex.axis=1.5,cex.main=2, main='Magnitude Error')
      graphics::axis(side= 1, at= plot.dates, labels= format(plot.dates, "%b %y"), cex.axis = 1.5)
      graphics::abline(h=0,lty=1, lwd=0.5, col='black')
      # end plots
      
      ########################################
      #  store segment-based stats in table  #
      ########################################
      #-------#
      #  obs  #
      #-------#
      #period
      obs.per = difftime(obs.segment$time[length(obs.segment$time)],obs.segment$time[1], units=time.units)
      
      #amplitude
      obs.amp           = obs.segment$val[length(obs.segment$val)] - obs.segment$val[1] 
      segment_type      = ifelse(obs.amp > 0, 'rise', 'fall')
     
      #segment start,end day of year
      obs.doy.start     = (as.POSIXlt(obs.segment$time[1], format = "%Y-%m-%d"))$yday
      obs.doy.end       = (as.POSIXlt(obs.segment$time[1], format = "%Y-%m-%d"))$yday

      #update stats
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'tracer']       <- obs.name
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'n_obs']        <- nrow(obs.segment)
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'segment_type'] <- segment_type
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'doy_start']    <- obs.doy.start
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'doy_end']      <- obs.doy.start
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'period']       <- abs(obs.per)
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'amplitude']    <- abs(obs.amp)
      
      #-------#
      #  sim  #
      #-------#
      #period
      sim.per = difftime(sim.segment$time[length(sim.segment$time)],sim.segment$time[1], units=time.units)
      
      #amplitude
      sim.amp = sim.segment$val[length(sim.segment$val)] - sim.segment$val[1] 
      segment_type = ifelse(obs.amp > 0, 'rise', 'fall')
      
      #segment start,end day of year
      sim.doy.start     = (as.POSIXlt(sim.segment$time[1], format = "%Y-%m-%d"))$yday
      sim.doy.end       = (as.POSIXlt(sim.segment$time[1], format = "%Y-%m-%d"))$yday

      #update stats
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'tracer']            <- sim.name
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'n_obs']             <- nrow(sim.segment)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'segment_type']      <- segment_type
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'doy_start']         <- sim.doy.start
      df_errorStats[which(df_errorStats$tracer==obs.name & df_errorStats$segment == seg_j),'doy_end']           <- sim.doy.end
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'period']            <- abs(sim.per)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'amplitude']         <- abs(sim.amp)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'timing_rmse']           <- sqrt(mean(segment.dist$dist_tdiff^2))
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'timing_bias_max']       <- max(segment.dist$dist_tdiff)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'timing_bias_min']       <- min(segment.dist$dist_tdiff)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'timing_bias_mean']      <- mean(segment.dist$dist_tdiff)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'timing_bias_var']       <- var(segment.dist$dist_tdiff)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'magnitude_rmse']      <- sqrt(mean(segment.dist$dist_vdiff^2))
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'magnitude_bias_max']  <- max(segment.dist$dist_vdiff)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'magnitude_bias_min']  <- min(segment.dist$dist_vdiff)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'magnitude_bias_mean'] <- mean(segment.dist$dist_vdiff)
      df_errorStats[which(df_errorStats$tracer==sim.name & df_errorStats$segment == seg_j),'magnitude_bias_var']  <- var(segment.dist$dist_vdiff)
    }#..end segment loop
    if(save.plot==TRUE){grDevices::dev.off()}
  #return data.frame for single region
  return(df_errorStats)
}#..end of fn

