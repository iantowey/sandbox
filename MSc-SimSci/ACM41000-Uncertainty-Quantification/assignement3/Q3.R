library(stringr)
library(dplyr)
library(R.utils)
library(R.cache)
library(Matrix)
library(fda)

#init env
rm(list=ls())
source("/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/Data2LD_Fix/Data2LD_Fix/R/make.coef.R")
source("/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/Data2LD_Fix/Data2LD_Fix/R/coefCheck.R")

#import data 
measles.data <- read.csv(file='/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/Measles.csv', header=TRUE, na.strings = c("*"), colClasses = rep("integer", 10))
num.years = measles.data %>% distinct(YY) %>% nrow
measles.data[is.na(measles.data)] <- 0

#sort data
measles.data = measles.data %>% mutate(dt=as.integer(paste(YY,str_pad(MM,2,pad="0"),str_pad(X.DD,2,pad="0"), sep=""))) %>% select(-c(1,2,3)) %>% arrange(-desc(dt)) %>% mutate(wk_num=row_number()) %>% select(-c(8))
measles.data <- measles.data[seq(1,316,1),]
measles.rng<- c(0,nrow(measles.data))
plot(measles.data$wk_num, measles.data$London, type = 'l')
plot(measles.data$wk_num, measles.data$Newcastle, type = 'l')

#set up vars
city.population <- list(London=8171000,Bristol=436000,Liverpool=747500, Manchester=661000, Newcastle=269400, Birmingham=1105000, Sheffield=494000)
conbasis <- create.constant.basis(measles.rng)

#create the coefficients
coefList  <- vector("list",length(names(city.population))^2)
i <- 1
for (city in names(city.population)){
  for (c in names(city.population)){
    coefList[[i]] <-make.coef(fun=conbasis, parvec=0.1, estimate=TRUE, coeftype="beta")
    i <- i + 1
  }
}
coefCheckList <- coefCheck(coefList)

# Define the terms in the homogeneous part of the equation
XListList = list()
i = 1
for (c in names(city.population)){
  j = 1
  XList <- vector("list",length(city.population))
  for (city in names(city.population)){
    XList[[j]] <- make.Xterm(variable=1, derivative=0, ncoef=i, factor=-1)
    i <- i + 1
    j <- j + 1
  }
  XListList[[c]] <- XList 
}

#  Define the system of differential equation in the model
measles.ODEList=vector("list",length(city.population))
idx = 1
for (city in names(city.population)){
  measles.ODEList[[idx]] <- make.variable(name=city,order=1, XList=XListList[[city]])
  idx <- idx + 1
}
#  check the object for internal consistency
measles.ODEList <- modelCheck(measles.ODEList, coefList)

#setup data
yList <- vector("list",7)
idx = 1
for (city in names(city.population)){
  yList[[idx]] <- list(argvals=seq(1,nrow(measles.data)), y=measles.data[[city]])
  idx <- idx + 1
}

#define basis functions
norder    <- 7
break.array <- seq(25,nrow(measles.data),52)
nbreaks   <- length(break.array)
nbasis    <- norder + nbreaks - 2
measles.Basis  <- create.bspline.basis(measles.rng, nbasis)
XbasisList <- vector("list",7)
XbasisList[[1]] <- measles.Basis
XbasisList[[2]] <- measles.Basis
XbasisList[[3]] <- measles.Basis
XbasisList[[4]] <- measles.Basis
XbasisList[[5]] <- measles.Basis
XbasisList[[6]] <- measles.Basis
XbasisList[[7]] <- measles.Basis

GMSLTimefine <- seq(0,nrow(measles.data),len=1001)
GMSLTBasisfine <- eval.basis(GMSLTimefine, measles.Basis)
matplot(GMSLTimefine, GMSLTBasisfine, type="l", xlim=c(0,nrow(measles.data)), ylim=c(0,1),
        xlab="Time t", ylab="Basis functions phi(t)")

dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg  <- c(1e-8, 1e-4)  #  convergence criterion
rhoVec <- rep(0.5,7)  #  light smoothing

length(yList[[1]]$y)
length(XbasisList) 
length(measles.ODEList) 
length(coefList) 
length(rhoVec)

Data2LDResult <- Data2LD(yList, XbasisList, measles.ODEList, coefList, rhoVec)
Data2LDResultOpt <- Data2LD.opt(yList, XbasisList, measles.ODEList, coefList,  rhoVec, convrg, iterlim, dbglev)
