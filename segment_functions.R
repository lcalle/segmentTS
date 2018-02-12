#============================================#
#            F U N C T I O N S               #
#     ::required for segment analysis::      #
#  wherenearest:                             #
#   ..nearest neighboor matching             #
#  fn_simplify_tracers:                      # 
#   ..main fn, makes plots, saves stats      #
#  SeriesDist.XX:                            #    
#   ..segment functions                      #
#============================================#
wherenearest        <- function(val,matrix){
  dist = abs(matrix-val)
  index = which.min(dist)
  return( index )
}
fn_simplify_tracers <- function(inDir,outDir,tracer_type, model,func.var,makePlot=FALSE,reg_j=reg_j, ...){

  #set up names
  region_plot = plot_title_codes[reg_j]
  model_name <- model

  #get data
  ccgcrv_g <- read.table(file = paste0(inDir,"mean_gosat_acos_for_ccgcrv_region_",region_plot,".txt"))
  ccgcrv_s <- read.table(file = paste0(inDir,"mean_",model_name,"_ACOS_for_ccgcrv_region_",region_plot,".txt"))
  ccgcrv_f <- read.table(file = paste0(inDir,"mean_fossilfuel_ACOS_for_ccgcrv_region_",region_plot,".txt"))
  ccgcrv_o <- read.table(file = paste0(inDir,"mean_ocn09_ACOS_for_ccgcrv_region_",region_plot,".txt"))
  
  # CURRENT NAMES for equal spaced returns in ccgcrv: '-equal' and omits output for '-orig' (actm_orig), '-detrend' (actm_origMINtrend), '-res' (fit_residuals), '-ressm' (fit_SmoothResiduals)                 
  names(ccgcrv_g)  <- names(ccgcrv_s) <- names(ccgcrv_f) <- names(ccgcrv_o) <- names_dfccgcrv
  
  #add dates
  ccgcrv_g$date_obj  <- as.Date(paste(ccgcrv_g$year,ccgcrv_g$month,ccgcrv_g$day,sep="-"),format="%Y-%m-%d")
  ccgcrv_s$date_obj  <- as.Date(paste(ccgcrv_s$year,ccgcrv_s$month,ccgcrv_s$day,sep="-"),format="%Y-%m-%d")
  ccgcrv_f$date_obj  <- as.Date(paste(ccgcrv_f$year,ccgcrv_f$month,ccgcrv_f$day,sep="-"),format="%Y-%m-%d")
  ccgcrv_o$date_obj  <- as.Date(paste(ccgcrv_o$year,ccgcrv_o$month,ccgcrv_o$day,sep="-"),format="%Y-%m-%d")
  
  #set event data and select start dates
  obs  <- ccgcrv_g[which((ccgcrv_g$date_obj %in% ccgcrv_g$date_obj) == TRUE),]
  sim  <- ccgcrv_s[which((ccgcrv_s$date_obj %in% obs$date_obj) == TRUE & (duplicated(ccgcrv_s$date_obj)==FALSE)),]
  
  #set empty data frame
  df.obs.evnt <- data.frame(val=obs[,func.var], time=obs$date_obj, ts=1, te=length(obs[,func.var]), match=0)
  df.sim.evnt <- data.frame(val=sim[,func.var], time=sim$date_obj, ts=1, te=length(sim[,func.var]), match=0)
  
  #--------------------------------
  # (1) match gosat and simulated
  #--------------------------------
  df.obs.evnt$match <- 1:df.obs.evnt$te[1]
  df.sim.evnt$match <- 1:df.sim.evnt$te[1]
  sim.ls.evnt       <- list(df.obs.evnt,df.sim.evnt)
  names(sim.ls.evnt) = c('obs.evnt', 'sim.evnt')      
  
  #--------------------------------
  # (2) determine curve positions {increase,decrease,peak,trough}
  #--------------------------------
  df.obs.valpos <- SeriesDist.2CurvePos(dat= df.obs.evnt, lolim= -9999)
  df.sim.valpos <- SeriesDist.2CurvePos(dat= df.sim.evnt, lolim= -9999)
  
  #--------------------------------
  # (3) equalize curve segments (remove false peaks,valleys)
  #--------------------------------
  ls.eq.evnt <- SeriesDist.3EqualizeEvents(obs.evnt= df.obs.valpos, sim.evnt= df.sim.valpos, region_plot=region_plot, is.dgvm= ifelse(tracer_type=="dgvm",model,'FALSE'))
  
  #----------------------------------------
  # Get distance data for single segment
  #----------------------------------------
  #select corresponding segment for OBS,SIMULATE
  obs.seg.time = ls.eq.evnt[['obs.eq']]
  sim.seg.time = ls.eq.evnt[['sim.eq']]
  
  #----------------------------------------
  # Get number of segments and Store Data
  #----------------------------------------
  num_segments_obs = nrow(obs.seg.time)-1
  num_segments_sim = nrow(sim.seg.time)-1

  #--------------------------------------------
  # loop over all tracers and get segment data
  #--------------------------------------------
  # ensure that there are segments for comparison
  # ..possible that more segments in obs than in sim, and vice versa
  for(foldcode_tracers in 1:1){
    num_segments = min(num_segments_obs, num_segments_sim)
    
    #empty data frame 
    df_errorStats= data.frame( region  = rep(plot_title_regions[reg_j], times=num_segments*2),
                               tracer  = rep(c('gosat',tracer_type),times=num_segments),
                               n_obs   = rep(0,num_segments*2),
                               segment = rep(1:num_segments, each=2),
                               timing_rmse      = rep(0,num_segments*2),
                               timing_bias_max  = rep(0,num_segments*2),
                               timing_bias_min  = rep(0,num_segments*2),
                               timing_bias_mean = rep(0,num_segments*2),
                               timing_bias_var  = rep(0,num_segments*2),
                               magnitude_rmse      = rep(0,num_segments*2),
                               magnitude_bias_max  = rep(0,num_segments*2),
                               magnitude_bias_min  = rep(0,num_segments*2),
                               magnitude_bias_mean = rep(0,num_segments*2),
                               magnitude_bias_var  = rep(0,num_segments*2),
                               period            = rep(0,num_segments*2),
                               amplitude         = rep(0,num_segments*2),
                               uptakeRelease_obs = rep(0,num_segments*2),
                               uptakeRelease_sim = rep(0,num_segments*2))
    
    if(makePlot==TRUE){
        pdf(paste0(outDir,"plot_smcycle_region_",plot_title_codes[reg_j],"_",model_name,"_seriesDistance.pdf"),
            height=10,width=18)
    }

    for(seg_j in 1:num_segments){
      seg_beg = seg_j
      seg_end = seg_j+1
      
      #----------------------------------
      # get event data for segment
      #----------------------------------
      obs.segment = df.obs.evnt[which(df.obs.evnt$time >= obs.seg.time$time[seg_beg] & df.obs.evnt$time <= obs.seg.time$time[seg_end]),]
      sim.segment = df.sim.evnt[which(df.sim.evnt$time >= sim.seg.time$time[seg_beg] & df.sim.evnt$time <= sim.seg.time$time[seg_end]),]
      
      #----------------------------------
      # get distance statistics
      #----------------------------------
      seriesDist.sim <- SeriesDist.4SegmentDistance(obs.seg = obs.segment , sim.seg = sim.segment)
      
      #############################
      #  P L O T : seriesDistance
      #############################
      #multi-plot with errors for segment
      #..only create one plot per model for the seasonal cycle
      if(makePlot==TRUE){
        par(mar=c(4,6,5,3),mfcol=c(2,2))

        #plot attributes
        x.lim = c(min(ccgcrv_g$date_obj,obs.segment$time[1],sim.segment$time[1]), max(obs.segment$time[length(obs.segment$time)],sim.segment$time[length(sim.segment$time)]))
        y.lim = c(min(ccgcrv_g[, func.var],obs.segment$val,sim.segment$val), max(ccgcrv_g[, func.var],obs.segment$val,sim.segment$val))
        plot.dates = seq.Date(from = x.lim[1], to = x.lim[2], by = 'month')
        plot.all.dates = seq.Date(from = ccgcrv_g[1,'date_obj'], to = ccgcrv_g[nrow(ccgcrv_g),'date_obj'], by = 'month')
        
        #plot
        plot(x= ccgcrv_g[,'date_obj'], y= ccgcrv_g[, func.var],lwd=3,
             type='l', col='black', ylab='XCO2 (ppm)', xlab='',
             xlim= c(x.lim[1], as.Date("2012-12-31")),
             ylim= y.lim, xaxt='n', cex.lab=2,cex.axis=1.5,
             cex.main=2,
             main=paste0('GOSAT (black);', ' DGVM ',model_name,' (red)\nXCO2 Seasonal Cycle in ',
             plot_title_regions[reg_j]))
        axis(side= 1, at= plot.all.dates, labels= format(plot.all.dates, "%b %y"), cex.axis = 1.5)
        lines(x= ccgcrv_s$date_obj,y= ccgcrv_s[,func.var], type='l', col='red', lwd=2)
        lines(x= ccgcrv_f$date_obj,y= ccgcrv_f[,func.var], type='l', col='grey50',lty=2, lwd=2)
        lines(x= ccgcrv_o$date_obj,y= ccgcrv_o[,func.var], type='l', col='grey50',lty=3, lwd=2)
      
        #points for matching segments
        points(x=obs.segment$time[seq(1,length(obs.segment$time), length=25)], y=obs.segment$val[seq(1,length(obs.segment$time), length=25)], pch=21,cex=1.2, bg='deepskyblue')
        points(x=sim.segment$time[seq(1,length(sim.segment$time), length=25)], y=sim.segment$val[seq(1,length(sim.segment$time), length=25)], pch=21,cex=1.2, bg='deepskyblue')

        abline(h=0,lty=1, lwd=0.5, col='black')

        #----------------------------------------
        # Plot Close-up of Matching Segments
        #----------------------------------------
        x.main       = ccgcrv_s[which((ccgcrv_s$date_obj %in% obs.segment$time)==TRUE),'date_obj']
        y.main       = ccgcrv_s[which((ccgcrv_s$date_obj %in% obs.segment$time)==TRUE), func.var]
        obs.segment  = obs.segment
        sim.segment  = sim.segment
        x.lim        = c(min(obs.segment$time[1],sim.segment$time[1]), max(obs.segment$time[length(obs.segment$time)],sim.segment$time[length(sim.segment$time)]))
        y.lim        = c(min(obs.segment$val,sim.segment$val), max(obs.segment$val,sim.segment$val))
        main.label   = paste0('Matching Segment')
        col.sim      = 'red'

        #close up plots
        plot(x= x.main, y= y.main,
             type='l', col=col.sim, ylab='ppm', xlab='',
             xlim= x.lim, ylim= y.lim, xaxt='n',
             cex.lab=2,cex.axis=1.5,
             cex.main=2, main=main.label)
        axis(side= 1, at= plot.dates, labels= format(plot.dates, "%b %y"), cex.axis = 1.5)
        lines(y= ccgcrv_g[,func.var], x= ccgcrv_g$date_obj,type='l')
        lines(x=obs.segment$time, y=obs.segment$val,type='l', col='black', lwd=4)
        lines(x=sim.segment$time, y=sim.segment$val,type='l', col=col.sim, lwd=4)
        points(x=obs.segment$time[seq(1,length(obs.segment$time), length=35)], y=obs.segment$val[seq(1,length(obs.segment$time), length=35)], pch=21,cex=1.2, bg='deepskyblue')
        points(x=sim.segment$time[seq(1,length(sim.segment$time), length=35)], y=sim.segment$val[seq(1,length(sim.segment$time), length=35)], pch=21,cex=1.2, bg='deepskyblue')
        #for clean figure, plot only some of the segments
        segments(x0= obs.segment$time[c(TRUE,rep(FALSE,5))], y0= obs.segment$val[c(TRUE,rep(FALSE,5))],
                 x1= seriesDist.sim[['poly_t']][c(TRUE,rep(FALSE,5))],
                 y1= seriesDist.sim[['poly_v']][c(TRUE,rep(FALSE,5))], col='grey75')
        abline(h=0,lty=1, lwd=0.5, col='black')

        #----------------------------------------
        # Plot Phase & Magnitude Error
        #----------------------------------------
        #plot phase error
        plot(x=obs.segment$time, xaxt='n', y=seriesDist.sim[['dist_tdiff']], type='o', pch=21, bg='grey25', ylab='Days', xlab='', xlim= x.lim,
             cex.lab=2, cex.axis=1.5, cex.main=2, main='Phase Eror')
        axis(side= 1, at= plot.dates, labels= format(plot.dates, "%b %y"), cex.axis = 1.5)
        abline(h=0,lty=1, lwd=0.5, col='black')
        #plot magnitude error
        plot(x=obs.segment$time, xaxt='n',  y=seriesDist.sim[['dist_vdiff']], type='o', pch=24, bg='grey25', ylab='ppm', xlab='', xlim= x.lim,
             cex.lab=2, cex.axis=1.5,cex.main=2,main='Magnitude Error')
        axis(side= 1, at= plot.dates, labels= format(plot.dates, "%b %y"), cex.axis = 1.5)       
        abline(h=0,lty=1, lwd=0.5, col='black')
        # end plots
      }#..end of foldcode
      
      #############################
      #  Store Data for Tracer
      #############################
      #-------------------------------------
      # calculate period of uptake/release
      # calculate amplitude of segment
      # ..store values for tracer type
      #-------------------------------------
      #----------------------------------
      # gosat - obs
      #----------------------------------
      #Carbon Uptake Period; CUP (days) [obs and sim]
      obs.cup = difftime(obs.segment$time[length(obs.segment$time)],obs.segment$time[1], units='days')
      
      #Amplitude (ppm) [Net: Source or Sink]
      obs.amp           = obs.segment$val[length(obs.segment$val)] - obs.segment$val[1] 
      uptakeRelease_obs = ifelse(obs.amp > 0, 'release', 'uptake')
      
      #update stats
      df_errorStats[which(df_errorStats$tracer=='gosat' & df_errorStats$segment == seg_j),'model']             <- 'gosat'
      df_errorStats[which(df_errorStats$tracer=="gosat" & df_errorStats$segment == seg_j),'uptakeRelease_obs'] <- uptakeRelease_obs
      df_errorStats[which(df_errorStats$tracer=="gosat" & df_errorStats$segment == seg_j),'period']            <- abs(obs.cup)
      df_errorStats[which(df_errorStats$tracer=="gosat" & df_errorStats$segment == seg_j),'amplitude']         <- abs(obs.amp)
      
      #----------------------------------
      # modeled - land,fossil, or ocean
      #----------------------------------
      #Carbon Uptake Period; CUP (days) [obs and sim]
      obs.cup = difftime(obs.segment$time[length(obs.segment$time)],obs.segment$time[1], units='days')
      sim.cup = difftime(sim.segment$time[length(sim.segment$time)],sim.segment$time[1], units='days')
      
      #Amplitude (ppm) [Net: Source or Sink]
      obs.amp = obs.segment$val[length(obs.segment$val)] - obs.segment$val[1] 
      sim.amp = sim.segment$val[length(sim.segment$val)] - sim.segment$val[1] 
      uptakeRelease_obs = ifelse(obs.amp > 0, 'release', 'uptake')
      uptakeRelease_sim = ifelse(obs.amp > 0, 'release', 'uptake')
      
      #update stats
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'model']             <- model_name
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'n_obs']             <- nrow(obs.segment)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'uptakeRelease_obs'] <- uptakeRelease_obs
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'uptakeRelease_sim'] <- uptakeRelease_sim
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'period']            <- abs(sim.cup)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'amplitude']         <- abs(sim.amp)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'timing_rmse']           <- sqrt(mean(seriesDist.sim$dist_tdiff^2))
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'timing_bias_max']       <- max(seriesDist.sim$dist_tdiff)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'timing_bias_min']       <- min(seriesDist.sim$dist_tdiff)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'timing_bias_mean']      <- mean(seriesDist.sim$dist_tdiff)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'timing_bias_var']       <- var(seriesDist.sim$dist_tdiff)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'magnitude_rmse']      <- sqrt(mean(seriesDist.sim$dist_vdiff^2))
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'magnitude_bias_max']  <- max(seriesDist.sim$dist_vdiff)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'magnitude_bias_min']  <- min(seriesDist.sim$dist_vdiff)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'magnitude_bias_mean'] <- mean(seriesDist.sim$dist_vdiff)
      df_errorStats[which(df_errorStats$tracer==tracer_type & df_errorStats$segment == seg_j),'magnitude_bias_var']  <- var(seriesDist.sim$dist_vdiff)
    }
    if(makePlot==TRUE){dev.off()}
    # end loop over segments ---------------------
  }#..end tracer loop

  return(df_errorStats)
}#..end of fn

#=======================================================================
# Functions for Series Distance (Ehret and Zehe 2011 HydroEarthSysSci)
# ..modified as in (Calle, Poulter, ana Patra 2018)
#=======================================================================
SeriesDist.1MatchEvents     <- function(obs.evnt = df, sim.evnt =  df, limit4match = 0){
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
SeriesDist.2CurvePos        <- function(dat = df, lolim = -999){  
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
SeriesDist.3EqualizeEvents  <- function(obs.evnt = df, sim.evnt = df,region_plot, is.dgvm="FALSE"){
  ##############################c
  #
  # D A T A   S T R U C T U R E
  # -- (obs.evnt, sim.evnt) data.frame with variables:
  # -- ($pkval) peak values, ordered by time
  # -- ($tgval) trough values, ordered by time
  #
  # Note: n= # of obs.evnts; m= # of sim.evnts
  # Returns: list of two data.frames, with number of peaks,troughs equalized
  # -- focuses on main signals, not local maxima/minima
  ##############################c
  #peak & trough values
  obs.peak   = obs.evnt[which(obs.evnt$pos== 2),]
  obs.trough = obs.evnt[which(obs.evnt$pos== -2),]
  sim.peak   = sim.evnt[which(sim.evnt$pos== 2),]
  sim.trough = sim.evnt[which(sim.evnt$pos== -2),]
  
  #-------------------------------------
  # match segment by visual inspection
  # ..rise-rise, decline-decline
  #-------------------------------------
  for(foldcode_manual_adjust_peakTrough in 1:1){
      if(region_plot==2){
        if(is.dgvm =="OCN" || is.dgvm =="VISIT"){
          obs.peak = obs.peak[c(-1),]
        }else if(is.dgvm =="LPJ"){
          obs.peak = obs.peak[c(-1),]
          sim.trough = sim.trough[c(-2,-4),]
        }else if(is.dgvm =="LPX"){
          obs.peak = obs.peak[c(-1),]
        }else if(is.dgvm =="JULES"){
          obs.peak   = obs.peak[c(-1),]
          sim.peak   = sim.peak[c(-1),]
          sim.trough = sim.trough[c(-1),]
        }
      }else if(region_plot==4){ 
        if(is.dgvm =="LPJ"){
          obs.peak = obs.peak[c(-1,-2,-3),]
          sim.peak = sim.peak[c(-3,-5),]
        }else if(is.dgvm =='JULES'){
          sim.trough = sim.trough[c(-1,-2),]
        }else if(is.dgvm =='OCN'){
          sim.peak   = sim.peak[c(-1),]
          sim.trough = sim.trough[c(-1),]
        }  
      }else if(region_plot==5){
        if(is.dgvm =="LPJ" || is.dgvm =="LPX" || is.dgvm =="VISIT"){
          obs.peak = obs.peak[c(-1),]
        }
      }else if(region_plot==6){
        if(is.dgvm == "CLM"){
          sim.peak = sim.peak[c(-1,-3),]
        }else if(is.dgvm == "ORCHIDEE"){
          sim.peak = sim.peak[c(-8),]
        }else if(is.dgvm == "OCN"){
          sim.peak = sim.peak[c(-6,-8),]
        }else if(is.dgvm == "LPX"){
          sim.peak = sim.peak[c(-3),]
        }else if(is.dgvm == "JULES"){
          sim.peak = sim.peak[c(-6,-8),]
          sim.trough = sim.trough[c(-1),]
        }
      }else if(region_plot==7){
        if(is.dgvm != "ORCHIDEE"){
          obs.peak = obs.peak[c(-1,-2),] 
        }
      }else if(region_plot==8){
        if(is.dgvm =="ORCHIDEE"){
          sim.peak = sim.peak[c(-1,-2),]
        }else if(is.dgvm =="CLM"){
          sim.peak = sim.peak[c(-1),]
        }
      }else if(region_plot==9){
        obs.peak = obs.peak[c(-2),]
        if(is.dgvm =="CLM" || is.dgvm =="ORCHIDEE"){
          sim.peak = sim.peak[c(-1),]
        }else if(is.dgvm =="JULES"){
          sim.peak   = sim.peak[c(-1,-3,-4),]
          sim.trough = sim.trough[c(-5),]
        }else if(is.dgvm =="OCN"){
          sim.peak = sim.peak[c(-1,-2,-3),]
        }
      }else if(region_plot==10){
        obs.peak = obs.peak[c(-1,-3,-5),]
        obs.trough = obs.trough[c(-1,-3,-5,-7),]
        if(is.dgvm =="JULES"){
          sim.peak = sim.peak[c(-1,-3),]
          sim.trough = sim.trough[c(-1,-3),] 
        }
      }else if(region_plot==11){
        if(is.dgvm =="LPJ" || is.dgvm =="LPX" || is.dgvm =="VISIT"){
          #sim.peak = sim.peak[c(-1),]
        }else if(is.dgvm=="CLM"){
          sim.peak = sim.peak[c(-1),]
        }
      }
  }#..end foldcode    

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
  #-------------------------------------------
  # try to equalize # of peaks, troughs 
  # ..manually adjust for regions, models
  #-------------------------------------------
  if(region_plot == 4 || region_plot == 6){
    val.mindays <- 150
  }else if(region_plot == 3 && is.dgvm=='JULES'){
    val.mindays <- 150
  }else{val.mindays <- 250}
  
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
      df = rbind(df.peak,df.trough)
      #sort by date/time
      df = df[order(df$time),]
      
      #keep only lowest of consective valleys and highest of consecutive peaks
      rm_vec = c()
      for(k in 1:(nrow(df)-1)){
        if(df$pos[k] == df$pos[k+1]){
          if(df$pos[k] == -2){
            #troughs: keep lowest value
            ifelse(df$val[k] < df$val[k+1], rm_vec <- c(rm_vec,-1*(k+1)), rm_vec <- c(rm_vec,-1*k))
          }else{
            #peaks: keep greatest value
            ifelse(df$val[k] > df$val[k+1], rm_vec <- c(rm_vec,-1*(k+1)), rm_vec <- c(rm_vec,-1*k))
          }
        }
      }
      #remove false peaks,valleys if they exist
      if(!is.null(rm_vec)){df = df[rm_vec,]}
      
      return(df)
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
    format_minDays_btwn_peakValley <- function(df, peakTrough = "peak or trough", min.days= 250){
      rm_val = NULL
      n= nrow(df)
      
      diff=1:(n-1)
      for(j in 1:(n-1)){
        diff[j]= abs(difftime(df$time[j],df$time[j+1], units='days'))
      } 
      
      j = which(diff < min.days,arr.ind = TRUE)
      if(length(j) > 0){
        rm_val=1:length(j)
        for(k in 1:length(j)){
          if(peakTrough == 'peak'){
            rm_val[k] = min(df$val[j[k]], df$val[j[k]+1])
          }else if(peakTrough == 'trough'){
            rm_val[k] = max(df$val[j[k]], df$val[j[k]+1])
          }
        }
      }#end do something if peaks/troughs closer than min.days
      
      #remove peak with lower value value
      if(!is.null(rm_val)){df <- df[!(df$val %in% rm_val),]}
      return(df)
}
#segment defined as peak to trough (rise or decline), equalizing peaks 'sim' by #peaks in 'obs'
#..stats are: simulation minus observation 
SeriesDist.4SegmentDistance <- function(obs.seg = df, sim.seg = df){
  ##############################c
  #
  # D A T A   S T R U C T U R E
  # -- takes in the segmented data (i.e., points within each unique $pos factor on curve)
  # -- (obs.seg, sim.seg) segement time-series data
  # -- ($val) values, ordered by time
  # -- ($time) timing of values on curve, ordered by time
  #
  # Note: n= # of obs.segs; m= # of sim.segs
  # Returns: two dfs; (1) dist.amp= (n) amplitude differences between obs and sims
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
