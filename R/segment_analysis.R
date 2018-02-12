#required functions
source("segmentTS/R/segment_functions.R")

#......................................................
#  M A I N   L O O P  -  S E R I E S  D I S T A N C E
#
# ..plot segment matching by model, region, segment 
# ..get statistics for series distance
#
# NOTE: ocean (ocn09) no real seasonal cycle,
#       ..so we omit below; code is place for ocn09
#......................................................
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
  saveRDS(object = list_errorStats_tracerRegion, file = paste0("segmentTS/outputs/seriesDistance_ccgcrv_region_",plot_title_codes[reg_j],"_",tracer_type,".RDS"))
}#..end region loop


