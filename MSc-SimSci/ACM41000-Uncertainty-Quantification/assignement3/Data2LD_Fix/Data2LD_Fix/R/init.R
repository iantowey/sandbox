init <- function(wd)

if(!is.loaded("RmatFn")){  
  paste(wd,"/src/RmatFn.so",sep="")
  dyn.load(paste(wd,"/src/RmatFn.so",sep=""))
}
if(!is.loaded("SmatFn")){ 
dyn.load(paste(wd,"/src/SmatFn.so",sep=""))
}
if(!is.loaded("loopJuan")){
dyn.load(paste(wd,"/src/loopJuan.so",sep=""))
}
if(!is.loaded("DRarrayFn")){
dyn.load(paste(wd,"/src/DRarrayFn.so",sep=""))
}
if(!is.loaded("DASarrayFn")){
dyn.load(paste(wd,"/src/DASarrayFn.so",sep=""))
}
if(!is.loaded("DBSarrayFn")){
dyn.load(paste(wd,"/src/DBSarrayFn.so",sep=""))
}

setwd(paste(wd,"/data",sep=""))
load(file = "HeadImpactData.rda")
load(file = "CanadianWeather.rda")
load(file = "fdascript.rda")
load(file = "RefineryData.rda")

setwd(paste(wd,"/R",sep=""))
