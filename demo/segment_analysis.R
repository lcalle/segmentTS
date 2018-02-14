#==================================================================
#  M A I N 
#
# OBJECTIVE: get segment-based statistics for two time-series,
#		   : ..create plots for visual inspection
# DATA     : use data(dfxco2). see man for details
#          : seasonal cycle time-series of XCO2
#          : ..for GOSAT (satellite) and Simulated XCO2
# FNS      : segment-based algorithms
#==================================================================
#---------------------------------
# define variables for analysis
#---------------------------------
#some region in the data
region.name = 'NA_Temperate'
#names of observation and simulated data
obs.name    = 'gosat'
sim.name    = 'CLM'

#name of variable for the data value we want in dfxco2 for seasonal cycle
func.var    = 'smcycle'

#name of date variable in dfxco2
date.var    = 'date_obj'

#pull sample data from package (redundant bc of lazyload)
data(dfxco2)

#-----------------------------#
#     P r e p  D a t a	      #
#-----------------------------#
#specify data and create the segment data frames 
#output is a list object w/ two data.frame objects
ls.evnt <- segmentTS.mkdf(df.obs=   dfxco2[which(dfxco2$model== obs.name & dfxco2$region_name== region.name),],
						  df.sim=   dfxco2[which(dfxco2$model== sim.name & dfxco2$region_name== region.name),],
						  func.var= func.var,
						  date.var= date.var)

#-----------------------------#
#         S t e p  1          #
#-----------------------------#
#specify the match index
#save as list object, updates data.frames from step1
#use default arg full.series=TRUE
ls.evnt <- segmentTS.1matchsignal(obs.evnt= ls.evnt[['obs.evnt']],
                                      sim.evnt= ls.evnt[['sim.evnt']],
                                      full.series=TRUE,limit4match = 0)

#-----------------------------#
#         S t e p  2          #
#-----------------------------#
#categorize signals
#returns positions of peak,trough,up,down,no-event
obs.evnt <- segmentTS.2catsignal(ls.evnt[['obs.evnt']], lolim = -999)
sim.evnt <- segmentTS.2catsignal(ls.evnt[['sim.evnt']], lolim = -999)

#-----------------------------#
#         S t e p  3          #
#-----------------------------#
#categorize signals
#returns list of two data.frames only with positions of equalized peaks,troughs
#after visual inspection of plots, if removal of peak/trough is neccessary, 
#..then specify the index of peak/trough, counting from left to right in plots
ls.evnt.pos <- segmentTS.3eqsignal.tmp(obs.evnt= obs.evnt,
                     sim.evnt= sim.evnt,
                     val.mindays = 250)

#-----------------------------#
#         S t e p  4          #
#-----------------------------#
#compute segment-based statistics and plots for all segments
#uses internal fn segmentTS.4segdist()
#returns a single data.frames with segment-based statistics
tmp.df <- segmentTS.4statsplots(obs.evnt= obs.evnt,
                      sim.evnt   = sim.evnt,
                      ls.evnt.pos= ls.evnt.pos,
                      obs.name   = obs.name,
                      sim.name   = sim.name,
                      time.units = 'days',
                      val.units  = 'ppm',
                      save.plot  = FALSE,
                      region.name= region.name)
