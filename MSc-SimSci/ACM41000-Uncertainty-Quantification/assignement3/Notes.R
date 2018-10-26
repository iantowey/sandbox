################################################
### initialise env
################################################
rm(list=ls())

library(R.utils)
library(R.cache)
library(Matrix)
library(fda)

wd <- "/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/Data2LD_Fix/Data2LD_Fix"
source("/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/Data2LD_Fix/Data2LD_Fix/R/init.R")
source("/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/Data2LD_Fix/Data2LD_Fix/R/make.coef.R")
source("/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/Data2LD_Fix/Data2LD_Fix/R/coefCheck.R")
init(wd)
v <- sourceDirectory("/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/Data2LD_Fix/Data2LD_Fix/R/")

##################################################
# Question 1
##################################################
data(ECBYieldCurve)
rate.ECB = first(ECBYieldCurve,'2 day')
maturity.ECB = c(0.25,0.5,seq(1,30,by=1))
NSParameters <- Nelson.Siegel(rate.ECB, maturity.ECB)
NS.rate <- NSrates(NSParameters,maturity.ECB)
plot(maturity.ECB, last(rate.ECB,'1 day'),main="Fitting Nelso Siegel yield curve",xlab=c("Pillars in years"), ylab=c("Rates"),type="o")
lines(maturity.ECB, last(NS.rate,'1 day'), col=2)
legend("topleft",legend=c("observed yield curve","fitted yield curve"),col=c(1,2),lty=1)
grid()


data(ECBYieldCurve)
rate.ECB = first(ECBYieldCurve,'2 day')
maturity.ECB = c(0.25,0.5,seq(1,30,by=1))
SvenssonParameters <- Svensson(rate.ECB, maturity.ECB)
Svensson.rate <- Srates( SvenssonParameters ,maturity.ECB,"Spot")
plot(maturity.ECB, last(rate.ECB,'1 day'),main="Fitting Svensson yield curve", xlab=c("Pillars in years"), ylab=c("Rates"),type="o")
lines(maturity.ECB, last(Svensson.rate,'1 day'), col=2)
legend("topleft",legend=c("observed yield curve","fitted yield curve"), col=c(1,2),lty=1)
grid()

##############
# (a)
##############

#define data
library(YieldCurve)
data(ECBYieldCurve)
rate.ECB <- as.vector(first(ECBYieldCurve,'1 day'))
length(rate.ECB)
maturity.ECB <- c(0.25,0.5,seq(1,30,by=1))
length(maturity.ECB)
YCrng <- c(0,max(maturity.ECB))
#define model
# Set up the constant basis
conbasis <- create.constant.basis(YCrng)
# Define the two constant coefficient functions
coef1 <- make.coef(fun=conbasis, parvec=0.25, estimate=TRUE, coeftype="beta")
coef2 <- make.coef(fun=conbasis, parvec=0.25,  estimate=TRUE, coeftype="beta")
#  Set up the coefficient list
coefList  <- vector("list",2)
coefList[[1]] <- coef1
coefList[[2]] <- coef2

#  check coefficient list
coefCheckList <- coefCheck(coefList)
coefList      <- coefCheckList$coefList
ntheta        <- coefCheckList$ntheta
print(paste("ntheta = ",ntheta))
# Define the two terms in the homogeneous part of the equation
Xterm0 <- make.Xterm(variable=1, derivative=1, ncoef=1, factor=-1)
Xterm1 <- make.Xterm(variable=1, derivative=2, ncoef=2, factor=-1)
# Set up the XList vector of length two
XList <- vector("list",2)
XList[[1]] <- Xterm0
XList[[2]] <- Xterm1

#  Define the single differential equation in the model
YCVariable <- make.variable(order=3, XList=XList)
YCVariableList=vector("list",1)
YCVariableList[[1]] <- YCVariable

#check model
YCVariableList <- modelCheck(YCVariableList, coefList)

#  Define the List array containing the data
yList1 = list(argvals=maturity.ECB, y=rate.ECB)
yList <- vector("list",1)
yList[[1]] <- yList1

plot(maturity.ECB, rate.ECB,main="yield curve",xlab=c("Pillars in years"), ylab=c("Rates"),type="o")

knots     <- 1:21
knots     <- c(0,0.5,0.5,0.5,1,1,1.5,1.5,2,2,3,3,4,5,seq(6,30,len=8))
length(knots)
norder    <- 6
nbasis    <- 26
YCBasis <- create.bspline.basis(YCrng,nbasis,norder,knots)

YCTimefine <- seq(0,30,len=201)
YCBasisfine <- eval.basis(YCTimefine, YCBasis)
matplot(YCTimefine, YCBasisfine, type="l", xlim=c(0,30), ylim=c(0,1), xlab="Time t", ylab="Basis functions phi(t)")

XbasisList <- vector("list",1)
XbasisList[[1]] <- YCBasis
# An evaluation of the criterion at the initial values
# Optimization of the criterion
#  algorithm constants
dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg  <- c(1e-8, 1e-4)  #  convergence criterion
rhoVec <- 0.995  #  light smoothing

Data2LDResult <- Data2LD(yList, XbasisList, YCVariableList, coefList, rhoVec)
Data2LDResultOpt <- Data2LD.opt(yList, XbasisList, YCVariableList, coefList,  rhoVec, convrg, iterlim, dbglev)

stderr  <- sqrt(diag( Data2LDResult$Var.theta)) 

Data2LDResultOpt$theta
thetaDn <- Data2LDResultOpt$theta - 2*stderr
thetaUp <- Data2LDResultOpt$theta + 2*stderr

cbind(Data2LDResultOpt$theta,thetaDn,thetaUp)

YCfd     <- Data2LDResult$XfdParList[[1]]$fd
YCfine   <- eval.fd(YCTimefine, YCfd)

plot(YCTimefine, YCfine, type="l",col="red", xlim=c(0,30),xlab="Time (msec)", ylab="Acceleration (cm/msec^2)")
points(maturity.ECB, rate.ECB, pch="o") 


Data2LDResultOpt$theta
MSE        <- Data2LDResult$MSE        # Mean squared error for fit to data
DpMSE      <- Data2LDResult$DpMSE      #  gradient with respect to parameter values
D2ppMSE    <- Data2LDResult$D2ppMSE    #  Hessian matrix
XfdParList <- Data2LDResult$XfdParList #  List of fdPar objects for variable values 
df         <- Data2LDResult$df         #  Degrees of freedom for fit
gcv        <- Data2LDResult$gcv        #  Generalized cross-validation coefficient
ISE        <- Data2LDResult$ISE        #  Size of second term, integrated squared error
SSE        <- MSE*df                   #  SSE
Var.theta  <- Data2LDResult$Var.theta  #  Estimate sampling variance for parameters

#  sequence of values of \rho
gammavec <- c(0:8)
rhoVec   <- exp(gammavec)/(1+exp(gammavec))
nrho     <- length(rhoVec)
#  Matrices to hold results
dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,1)
ISEsave <- matrix(0,nrho,1)
thesave <- matrix(0,nrho,ntheta)
#  Initialize coefficient list
coefList.opt <- coefList

for (irho in 1:nrho) {
  rhoi <- rhoVec[irho]
  print(paste(" ------------------  rhoVeci <- ", round(rhoi,4), 
              " ------------------"))
  Data2LDResult <- Data2LD.opt(yList, XbasisList, YCVariableList, coefList.opt, 
                               rhoi, convrg, iterlim, dbglev)
  theta.opti     <- Data2LDResult$theta
  coefList.opti  <- modelVec2List(theta.opti, coefList)
  DataListResult <- Data2LD(yList, XbasisList, YCVariableList, coefList, rhoVec)
  thesave[irho,] <- theta.opti
  dfesave[irho]  <- DataListResult$df
  gcvsave[irho]  <- DataListResult$gcv
  MSEsave[irho]  <- DataListResult$MSE
  ISEsave[irho]  <- DataListResult$ISE
  coefList.opt   <- coefList.opti
}




samp <- 50
d <- 3
mus <- rbind(rep(0,d), rep(1,d))
Sigmas <- 0.25*rbind(diag(d), diag(d))
props <- c(2/3, 1/3)
x <- rmvnorm.mixt(samp, mus, Sigmas, props)
H <- Hpi(x)
ise.mixt(x, H, mus, Sigmas, props)
mise.mixt(H, mus, Sigmas, props, samp)