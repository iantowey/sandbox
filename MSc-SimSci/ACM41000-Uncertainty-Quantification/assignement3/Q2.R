data <- read.csv(file='/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignement3/GMSL.csv', header=TRUE)
dim(data)
head(data, n=30)
plot(data, type = 'l')
rng<- c(0,nrow(data))

data$idx <- rownames(data)
GMSL.min <- -min(data$GMSL)
data = data %>% mutate(year = floor(Date), rescale.GMSL = GMSL + GMSL.min)
data.max.tmp.by.year =  data %>% group_by(year) %>% summarise(max.rescale.GMSL = max(rescale.GMSL))

pulse.idx = c(0,as.integer((data %>% inner_join(data.max.tmp.by.year, by = "year") %>% filter(rescale.GMSL == max.rescale.GMSL) %>%select(idx))$idx))


which.max(data[data$Date < 1994,]$GMSL)
which.max(data[data$Date >= 1994 & data$Date < 1995 ,]$GMSL)

conbasis <- create.constant.basis(HeadImpactRng)
# Define the three constant coefficient functions
coef1 <- make.coef(fun=conbasis, parvec=0.1, estimate=TRUE, coeftype="beta")
coef2 <- make.coef(fun=conbasis, parvec=0.1,  estimate=TRUE, coeftype="alpha")

coefList  <- vector("list",2)
coefList[[1]] <- coef1
coefList[[2]] <- coef2

coefCheckList <- coefCheck(coefList)
coefList      <- coefCheckList$coefList
ntheta        <- coefCheckList$ntheta
print(paste("ntheta = ",ntheta))
# Define the two terms in the homogeneous part of the equation
Xterm0 <- make.Xterm(variable=1, derivative=0, ncoef=1, factor=-1)
# Set up the XList vector of length two
XList <- vector("list",1)
XList[[1]] <- Xterm0

Pulsebasis <- create.bspline.basis(rng, 27, 1, pulse.idx)
Pulsefd    <- fd(matrix(c(0,1,0),27,1),Pulsebasis)

Fterm      <- make.Fterm(ncoef=2, Ufd=Pulsefd, factor=1)
# Set up the forcing term list of length one
FList      <- vector("list",1) 
FList[[1]] <- Fterm

#  Define the single differential equation in the model
GMSL.ODE <- make.variable(order=1, XList=XList, FList=FList)

GMSL.ODEList=vector("list",1)
GMSL.ODEList[[1]] <- GMSL.ODE
#  check the object for internal consistency
GMSL.ODEList <- modelCheck(GMSL.ODEList, coefList)

yList1 = list(argvals=data$Date, y=data$GMSL)
yList <- vector("list",1)
yList[[1]] <- yList1

knots     <- 0:nrow(data)
norder    <- 6
nbasis    <- 899
GMSL.Basis <- create.bspline.basis(rng,nbasis,norder,knots)

GMSLTimefine <- seq(0,nrow(data),len=1001)
GMSLTBasisfine <- eval.basis(GMSLTimefine, GMSL.Basis)
matplot(GMSLTimefine, GMSLTBasisfine, type="l", xlim=c(0,nrow(data)), ylim=c(0,1),
        xlab="Time t", ylab="Basis functions phi(t)")
lines(c(14,14), c(0,1), lty=2) 
lines(c(15,15), c(0,1), lty=2)

XbasisList <- vector("list",1)
XbasisList[[1]] <- GMSL.Basis


dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg  <- c(1e-8, 1e-4)  #  convergence criterion
rhoVec <- 0.995  #  light smoothing

Data2LDResult <- Data2LD(yList, XbasisList, GMSL.ODEList, coefList, rhoVec)
Data2LDResultOpt <- Data2LD.opt(yList, XbasisList, GMSL.ODEList, coefList,  rhoVec, convrg, iterlim, dbglev)
