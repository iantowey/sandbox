###  Analyses of motorcycle impact data

#  Last modified 11 November 2015

#  add path to fda functions

library(fda)

#  bypass set-up phase if previously completed

#  Two analyses: 
#    -- Using the full data and a pulse forcing function (here)
#    -- Using the post-impact data and no forcing function  (line 400)

##               Post- and Pre- impact Analysis


## Set up the full data

#  load the data

temp = scan("motorcycledata.txt",0)
motorcycledata = matrix(temp,133,3,byrow=TRUE)

motot = motorcycledata[,2]  #  time in milliseconds
motoy = motorcycledata[,3]  #  deformation in 0.1 millimeters

#  adjust the data for baseline and plot

impact  = 13.7  #  impact time
baseind = motot < impact
basey   = mean(motoy[baseind])

#  remove baseline, change time, and convert to centimeters

motoy   = (basey - motoy)/100.0  

yListi = list(argvals=motot, y=motoy)

yList = vector("list",1)
yList[[1]] = yListi

## plot data

plot(motot,motoy, type="p", xlim=c(0,60), ylim=c(-1.0,1.5),
     xlab="Time (milliseconds)", ylab="Acceleration (cm/msec^2)")

## Set up the basis system
#  Order 6 spline basis, with three knots at the impact point and 
#  three knots at impact + delta to permit discontinuous first 
#  derivatives at these points

motorng   = c(0,60)
knots = c(0,impact,impact,impact,impact+1,impact+1,seq(impact+1,60,len=11))
norder    = 6
delta     = 1
nbasis    = 21
motobasis = create.bspline.basis(motorng,nbasis,norder,knots)

basisList = vector("list",1)
basisList[[1]] = motobasis

## Set up the List array coefList

Wbasisobj = create.constant.basis(motorng)
WfdParobj = fdPar(Wbasisobj, 0, 0, TRUE)

coefList  = vector("list", 3)
coefList1 = list(fdPar = WfdParobj, index = 1)
coefList2 = list(fdPar = WfdParobj, index = 2)
coefList3 = list(fdPar = WfdParobj, index = 3)
coefList[[1]] = coefList1
coefList[[2]] = coefList2
coefList[[3]] = coefList3

coefcheckList = coefcheck(coefList)
coefList = coefcheckList$coefList
ntheta    = coefcheckList$ntheta

print(paste('ntheta = ',ntheta))

## Set up the List array modelList

Xterm0List = list(variable = 1, derivative = 0, ncoef = 2)
Xterm1List = list(variable = 1, derivative = 1, ncoef = 2)
XList = vector("list",2)
XList[[1]] = Xterm0List
XList[[2]] = Xterm1List

Ubasis = create.bspline.basis(motorng, 3, 1, c(0,impact,impact+1,60))
Ufd = fd(matrix(c(0,1,0),3,1),Ubasis)
Fterm = list(ncoef = 3, Ufd = Ufd)
FList = vector("list",1)
FList[[1]] = Fterm

order = 2
nallXterm = 2
nallFterm = 1
modelListi = list(XList = XList, FList = FList, order = order, 
                  nallXterm = nallXterm, nallFterm = nallFterm)

#  modelListi = list(XList = XList, FList = FList, order = order)

modelList = vector("list", 1)
modelList[[1]] = modelListi

modelList = modelcheck(modelList, coefList)

## An evaluation of the criterion at the initial values

rhoVec = 0.5  #  light smoothing
print(paste('rhoVec and 1 - rhoVec:  ',c(1-rhoVec,rhoVec)))
print(paste('rhoVec/(1-rhoVec) =     ',rhoVec/(1-rhoVec))))

[SSE, DSSE, D2SSE, XfdParList.opt, df, gcv, ISE] = 
                     D2LD(yList, basisList, 
                              modelList, coefList, rhoVec)

motofd.0 = getfd(XfdParList.opt[[1]])

## Optimize the fit for a range of rhoVec-values

rhoVec     = [0.50, 0.73, 0.88, 0.95, 0.98, 0.99, 0.999, 0.9999]
nrho       = length(rhoVec)
ntheta     = length(theta.opt)
thetastore = zeros(ntheta,nrho)
dfstore    = zeros(ntheta,1)
gcvstore   = zeros(ntheta,1)

for irho = 1:nrho
    rhoVeci = rhoVec(irho)
    theta.opti = D2LD.Opt(yList, basisList, modelList, coefList, 
                             rhoVeci, conv, iterlim, dbglev)
    coefList.opti = BAwtvec2List(theta.opti, coefList)
    [SSE, DSSE, D2SSE, XfdParList, df, gcv, ISE] = 
                   D2LD(yList, basisList, modelList, coefList.opti, rhoVec)
    thetastore(:,irho) = theta.opti
    dfstore(irho)      = df
    gcvstore(irho)     = gcv
end

# display the optimal parameter values

print(paste('Stiffness = ',thetastore(1,:))) 
print(paste('Damping   = ',thetastore(2,:)))
print(paste('Forcing   = ',thetastore(3,:)))

# Stiffness = -0.076   -0.077   -0.077   -0.073   -0.067   -0.063
# Damping   =  0.131    0.102    0.063    0.024   -0.029   -0.086
# Forcing   =  0.300    0.278    0.213    0.116   -0.032   -0.215
 
# period = 9.3 msec for 4th stiffness coef.

# display degrees of freedom and gcv values

print(paste(rhoVec',dfstore, gcvstore])

## Evaluate the fit for the minimum-gcv parameter values

irho = 6
theta4 = thetastore(:,irho)
coefList4 = BAwtvec2List(theta4, coefList)
rhoVec4 = rhoVec(irho)
[SSE, DSSE, D2SSE, XfdParList, df, gcv, ISE, Rvar, 
          coef, Dcoef, Cmat, Rmat, Smat] = 
                     D2LD(yList, basisList, modelList, coefList4, rhoVec4)

print(paste('SSE = ', SSE))
print(paste('df  = ', df))
print(paste('gcv = ', gcv))
print(paste('DSSE = ',DSSE'))

motofd.4 = getfd(XfdParList[[1]])

tfine = [linspace(motorng(1),impact,2),  
         linspace(impact,impact+delta,11), 
         linspace(impact+delta,motorng(2),101)]'
     
motovec.0 = eval.fd(tfine, motofd.0)
motovec.4 = eval.fd(tfine, motofd.4)

figure(3)
phdl=plot(tfine, motovec.0, 'r-', tfine, motovec.4, 'b-')
set(phdl, 'LineWidth', 2)
legend('\fontsize{16} \rho = 0.50', '\fontsize{16} \rho = 0.95', 
       'Location', 'SouthWest')
hold on
phdl=plot(motot, motoy, 'bo', 
          [impact,      impact],       [0,1], 'b--', 
          [impact+delta,impact+delta], [0,1], 'b--', 
          [impact,      impact+delta], [1,1], 'b--', 
          [0,60], [0,0], 'b:')
set(phdl, 'LineWidth', 2)
hold off
axis([0,60,-1.0,1.5])
xlabel('\fontsize{16} Time (milliseconds)')
ylabel('\fontsize{16} Acceleration (cm/msec^2)')
text(40,1.3,'\fontsize{16} \rho   =  0.95')
text(40,1.1,'\fontsize{16} \beta.0 = -0.073')
text(40,0.9,'\fontsize{16} \beta.1 =  0.024')
text(40,0.7,'\fontsize{16} \alpha  =  0.116')
text(40,0.5,'\fontsize{16} df  =  11.7')