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
inDir = "data/"
outDir= "outputs/"
#set region codes and corresponding region names for reference
plot_title_codes   <- c(2,3,4,11,5,6,8,7,9,10)
plot_title_regions <- c("NA_Temperate","SA_Tropical", "SA_Temperate", "Europe","Africa_Northern","Africa_Southern","Eurasia_Temperate","Eurasia_Boreal","Asia_Tropical","Australia")
tracers = c('dgvm'); dgvms = c("CLM", "JULES", "LPJ", "LPX", "OCN", "ORCHIDEE","VISIT") 

for(reg_j in 1:length(plot_title_codes)){
  #initlize list, counter
  list_errorStats_tracerRegion = list()
  count = 1
  
  #get data 
  for(dgvm in dgvms){
    #-------------------------------------------------
    # STORE stats for gosat, dgvm, fossilfuel
    #-------------------------------------------------
    list_errorStats_tracerRegion[[count]] <- fn_simplify_tracers(inDir=inDir,outDir=outDir,
                                                                 tracer_type = tracer_type,
                                                                 model = dgvm,
                                                                 func.var = "smcycle",
                                                                 makePlot=TRUE,
                                                                 reg_j=reg_j)
    #update counter
    count = count+1
  }

  #save data bby region
  saveRDS(object = list_errorStats_tracerRegion, file = paste0(outDir,"seriesDistance_ccgcrv_region_",plot_title_codes[reg_j],"_",tracer_type,".RDS"))
}#..end region loop


