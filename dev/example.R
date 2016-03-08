##Sample R code to testing SPRATS package
###By Feng 02/23/2016

#########testing code#######
dt<-new("SensorgramData",
	#associationData=data.frame(time=1:5,RU=1:5), 
	dissociationData=data.frame(time=1:5,RU=1:5, time2=1:5, RU=1:5))
dt@associationData
setwd("E:\\feng\\LAB\\hg\\SPR_models\\SPRTwoStatePackage");
dt<-read.table("CBSInhibitor_140313Run_5E-5ConcChanel1_raw.txt", skip=4, 
header=T, sep="\t")
associationPhaseStart<-0
associationPhaseEnd<-90
dissociationPhaseEnd<-120
x_rawData<-ReadSensorgramData("CBSInhibitor_140313Run_5E-5ConcChanel1_raw.txt", skip=4, 
header=T, sep="\t", associationPhaseEnd=80, dissociationPhaseEnd=150)

lgm<-new("LangmuirModel", kon=2E3, koff=0.001, analyteConcentrations=c(1E-6, 2E-6, 1E-5), associationLength=100, 
	dissociationLength=100,Rmax=50)
lgm2<-LangmuirModel(1e3, 0.01, c(1e-6, 5e-6, 1e-5), 30, 1000, 1000)	
lgIF<-new("InducedFitModel", BaseModel=lgm, kr=0.1, kf=0.1)
lgIF2<-InducedFitModel(1e3, koff=0.01,kr=0.01,kf= 0.01, analyteConcentrations=c(1e-6, 5e-6, 1e-5), 
	Rmax=50, associationLength=1000, dissociationLength=1000)	
	
lgCS<-new("ConformationalSelectionModel", BaseModel=lgm, kr=0.1, kf=0.1)
lgCS2<-ConformationalSelectionModel(1e3, koff=0.01,kr=0.01,kf= 0.01, analyteConcentrations=c(1e-6, 5e-6, 1e-5), 
	Rmax=100, associationLength=1000, dissociationLength=1000)	

#for two state model
lgTS<-new("TwoStateModel", BaseModel1=lgIF, BaseModel2=lgCS)	

lgTS2<-TwoStateModel( 1e5, koff=0.005,kr=0.002,kf= 0.003,
	kon2=1e4, koff2=0.001,kf2=0.01,kr2= 0.002,
	analyteConcentrations=c(5e-7,1e-6, 2e-6, 4e-6,8e-6), 
	Rmax=250, associationLength=1000, dissociationLength=1000)
	
 xt<-Simulate(lgTS2, timeStep=0.002, sampleFreq=0.5)

dt_AB<-xt[[1]]
dt_AB_star<-xt[[2]]	

#add the things together association + dissociation data
dt_RUs<-GetObservedRUs(dt_AB, dt_AB_star)
plot(dt_RUs)

kA1<-0.005/1E5
kA2<-0.001/1E4

kD_m<-kA1*0.002/(0.01+0.002)+kA2*0.01/(0.01+0.002)

koff<-(0.001*0.003+0.005*0.002)/0.005

#############new code section
	#first read in the data
setwd("E:\\feng\\LAB\\hg\\SPR_models\\SPRTwoStatePackage");
setwd("E:\\feng\\LAB\\MSI\\SPR\\SPRTwoStateData");
	dt_rd<-ReadSensorgramData("1558-steadystate_export.txt", skip=4, 
header=T, sep="\t", associationPhaseEnd=2790, dissociationPhaseEnd=3790)

	dt_rd<-ReadSensorgramData("Langmuir_wrong_sila.txt", skip=4, 
header=T, sep="\t", associationPhaseEnd=1000, dissociationPhaseEnd=2000)

	#dt_rd@analyteConcentrations<-c(9.60E-08,1.44E-07, 3.24E-07,2.16E-07,4.88E-07)
	dt_rd@analyteConcentrations<-c(1e-06,	2e-06,	1e-05,	2e-05,	1e-04)
	dt_rd@analyteConcentrations<-c(0.00E+00,9.60E-08,1.44E-07,2.17E-07,	0.00E+00,3.25E-07,4.88E-07,0.00E+00)

	dt_rd@steadyStateStart<-2700
	dt_rd@steadyStateEnd<-2790
	
	ret<-FitTwoStateSPR(dt_rd,type=1, init.association=list(Rmax=250,KD=5E-5),trace=T)
	ret_MeanSPR<-ret$MeanSPR
	ret_Concs<-ret$Concs
	fit_ret<-nls(ret_MeanSPR~Rmax*ret_Concs/(KD+ret_Concs), start=list(Rmax=250, KD=5E-7))
	rev_MeanSPR<-1/ret$MeanSPR
	rev_Concs<-1/ret$Concs
	fit_lm<-lm(rev_MeanSPR~rev_Concs)
	conc_by_Rss<-ret$Concs/ret$MeanSPR
	fit_lm_2<-lm(ret$Concs~conc_by_Rss)
	
	plot(c(1,100), c(90,150),type="n")
	for(i in c(1:8))
	{
		lines(dt_rd@dissociationData[c(1,100),i*2-1],dt_rd@dissociationData[c(1,100),i*2], type="l")
	}
	x<-dt_rd
	i<-2
	windowSize<-10
	times<-x@dissociationData[,i*2-1]
					RUs<-x@dissociationData[,i*2]
					tms<-times[times<windowSize]
					rus<-RUs[times<windowSize]