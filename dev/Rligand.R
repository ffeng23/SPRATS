##this is the R code to testing variable ligand levels model.

##first try to do simulation

setwd("E:\\feng\\LAB\\hg\\SPR_models\\SPRATS\\dev")
#Langmuir Model
lgm<-new("LangmuirModel", kon=2E3, koff=0.001, analyteConcentrations=c(1E-6, 2E-6, 1E-5), associationLength=200, 
	dissociationLength=100,Rmax=50)
lgm@Rligand<-c(50,50,50)
lgm2<-new("LangmuirModel", kon=2E3, koff=0.001, analyteConcentrations=c(2E-6, 4E-6, 8E-5, 1.5E-5), associationLength=200, 
	dissociationLength=1000,Rmax=60)
lgm2@Rligand<-c(60,50,40,30)
x<-Simulate(lgm, sampleFreq=1)

x2<-Simulate(lgm, sampleFreq=1, fix.ligand=FALSE)
plot(x[[1]])


#induced fit
lgIF<-new("InducedFitModel", BaseModel=lgm2, kr=0.1, kf=0.1)

xIF<-Simulate(lgIF, sampleFreq=1)
plot(xIF[[1]])
xIF2<-Simulate(lgIF, sampleFreq=1, fix.ligand=FALSE)

#CS
lgCS<-new("ConformationalSelectionModel", BaseModel=lgm2, kr=0.1, kf=0.1)
xCS<-Simulate(lgCS, sampleFreq=1)
plot(xCS[[1]])
xCS2<-Simulate(lgCS, sampleFreq=1, fix.ligand=FALSE)

#twostate
lgTS<-new("TwoStateModel", BaseModel1=lgIF, BaseModel2=lgCS)	
xts<-Simulate(lgTS, timeStep=0.002, sampleFreq=0.5)
plot(xts[[1]])
xts2<-Simulate(lgTS, timeStep=0.002, sampleFreq=0.5, fix.ligand=FALSE)
plot(xts2[[1]])

#####ready to do the fitting
dt_AB<-xts[[1]]
dt_AB_star<-xts[[2]]	

#add the things together association + dissociation data
dt_RUs<-GetObservedRUs(dt_AB, dt_AB_star)

plot(dt_RUs)

save(dt_RUs, file="simulatedTwoStateVar.RData")
load("simulatedTwoStateVar.RData");
#calling the fitting on variable ligand model, but testing the fix ligand model.
#of course this should work 
ret<-FitTwoStateSPR(dt_RUs,type=2, steadyStateStart=1800,steadyStateEnd=2000, windowSize=100,
					init.association=list(efficiency=1.2,KD=1E-5),trace=T, Rligand=c(60,60,60,60), fix.ligand=FALSE )


dt_AB2<-xts2[[1]]
dt_AB_star2<-xts2[[2]]	

#add the things together association + dissociation data
dt_RUs2<-GetObservedRUs(dt_AB2, dt_AB_star2)
plot(dt_RUs2)
save(dt_RUs2, file="simulatedTwoStateVar_variableLigand.RData")
load("simulatedTwoStateVar_variableLigand.RData");
#now testing the variable ligand model on variable ligand simulated data.
ret2<-FitTwoStateSPR(dt_RUs2,type=2, steadyStateStart=1800,steadyStateEnd=2000, windowSize=100,
					init.association=list(efficiency=1.2,KD=1E-5),trace=T, Rligand=c(60,50,40,30), fix.ligand=FALSE )
#write the data to 
xt<-dt_RUs2@dissociationData
tm4<-xt[,7]
dt4<-xt[,8]
dt4<-dt4+rnorm(length(tm4),0,0.1)
nlr2<-nls(dt4~R1*exp(r1*tm4)+R2*exp(r2*tm4), start=list(R1=-0.023, R2=29,r1=-0.0078,r2=-0.001),
							control=list(maxiter = 700,tol = 1e-2, minFactor=1/1E10), 
							trace = TRUE)
