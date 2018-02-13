#============================================================================================
# Fns for segmentTS as in (Calle, Poulter, ana Patra 2018)
#============================================================================================

#' Make data frames for event/signal type of time-series data
#'
#' This function makes the empty data frames for storing time-series event/signal data.
#' Removes data in simulated dataset that does not overlap (in time) with observed data.
#' @param df.obs data.frame object based on ccgcrv output for observational data.
#' At minimum, requires variables of values and date (YYYY-MM-DD)
#' @param df.sim data.frame object based on ccgcrv output for simulated data
#' At minimum, requires variables of values and date (YYYY-MM-DD)
#' @param func.var name of variable for data values. Must have same name in both df.obs and df.sim.
#' @param date.var name of variable for dates. Must have same name in both df.obs and df.sim.
#' @return a list object with two data.frames, df.obs.evnt, df.sim.evnt formatted with new variable names
#' @export
segmentTS.mkdf <- function(df.obs,df.sim,func.var,date.var){

  #make sure there is overlap in dates between sim and obs datasets
  df.sim  	  <- df.sim[which((df.sim[,date.var] %in% df.obs[,date.var]) == TRUE & (duplicated(df.sim[,date.var])==FALSE)),]

  #set empty data frame
  df.obs.evnt <- data.frame(val=df.obs[,func.var], time=df.obs[,date.var], ts=1, te=length(df.obs[,func.var]), match=0)
  df.sim.evnt <- data.frame(val=df.sim[,func.var], time=df.sim[,date.var], ts=1, te=length(df.sim[,func.var]), match=0)

  return(list(df.obs.evnt=df.obs.evnt,df.sim.evnt=df.sim.evnt))
}
