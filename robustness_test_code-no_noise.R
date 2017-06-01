library(ADASPR)
kon<-c(1E3, 1E3, 1E4)
koff<-c(1E-3, 5E-4, 1E-5)
analyteConcentrations<-c(3E-5, 1.5E-5, 1E-5, 7E-6, 5E-6)
associationLength<-1500 
#associationLength<-1000 
#associationLength<-500 
#associationLength<-100 
dissociationLength<-1500
Rmax<-c(100, 80, 60)


mlgm<- new("MultiLigandModel", kon=kon, koff=koff, analyteConcentrations=analyteConcentrations, 
           associationLength=associationLength, dissociationLength=dissociationLength, Rmax=Rmax)
set.seed(2)		
sData<-Simulate(mlgm,sampleFreq=0.1, sd=0)	 #for kon us sampleFreq=0.05
plot(sData[[1]])
fss<-FitSteadyStateSPR(sData[[1]], degree=5, steadyStateStart=1490,steadyStateEnd=1500, auto=T)
#fss<-FitSteadyStateSPR(sData[[1]], degree=5, steadyStateStart=990,steadyStateEnd=1000, auto=T)
#fss<-FitSteadyStateSPR(sData[[1]], degree=5, steadyStateStart=490,steadyStateEnd=500, auto=T)
#fss<-FitSteadyStateSPR(sData[[1]], degree=5, steadyStateStart=90,steadyStateEnd=100, auto=T)

fpc.on<-fitSPR.kon(sData[[1]],debug=TRUE,weights.type="exp", degree=7, weights.step=0.05, weights.scale=83,#25
			Rmax=240) #,mode=2
			#step 10~20, weights.scale=8.5


e_k<-rep(0,length(analyteConcentrations))
e_k[1]<-sum(koff/kon*(Rmax/sum(Rmax)))
e_k[2]<-sum((koff/kon)^2*(Rmax/sum(Rmax)))
e_k[3]<-sum((koff/kon)^3*(Rmax/sum(Rmax)))
e_k[4]<-sum((koff/kon)^4*(Rmax/sum(Rmax)))

#first get the distribution of Rmax
E_kon<-rep(0,length(kon))
E_kon[1]<-sum(Rmax/sum(Rmax)*kon)
E_kon[2]<-sum(Rmax/sum(Rmax)*kon^2)
E_kon[3]<-sum(Rmax/sum(Rmax)*kon^3)
E_kon[4]<-sum(Rmax/sum(Rmax)*kon^4)
E_kon[5]<-sum(Rmax/sum(Rmax)*kon^5)

fpc.off<-fitSPR.koff(sData[[1]], debug=TRUE,degree.fitMoments=4,  weightsType.fitSPR="uniform",degree.fitSPR=6, weightsScale.fitSPR=220#220, 
			,weightsStep.fitSPR=0.5#0.05, 
			,weightsType.fitMoments="exp", weightsScale.fitMoments=0.5);
			
for(i in 1:length(analyteConcentrations))
{
	E_koff1[i]<-sum(koff*r0[i,]/sum(r0[i,]))
	E_koff2[i]<-sum(koff^2*r0[i,]/sum(r0[i,]))
	E_koff3[i]<-sum(koff^3*r0[i,]/sum(r0[i,]))
	E_koff4[i]<-sum(koff^4*r0[i,]/sum(r0[i,]))
	
}


#expected koff moments over Rmax
E_koff<-rep(0,7)
E_koff[1]<-sum(Rmax)
E_koff[2]<-sum(koff*Rmax/sum(Rmax))
E_koff[3]<-sum(koff^2*Rmax/sum(Rmax))
E_koff[4]<-sum(koff^3*Rmax/sum(Rmax))
E_koff[5]<-sum(koff^4*Rmax/sum(Rmax))
E_koff[6]<-sum(koff^5*Rmax/sum(Rmax))
E_koff[7]<-sum(koff^6*Rmax/sum(Rmax))
			
