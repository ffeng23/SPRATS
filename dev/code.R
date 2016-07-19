lm<-LangmuirModel(kon=1E5, koff=1E-3, analyteConcentrations=c(1E-9, 2E-9, 1E-8), sd=0.5, Rmax=250, associationLength=500, dissociationLength=500)
ifm<-InducedFitModel(kon=1.72e5, koff=12.5e-4,kr=1.34e-4,kf=40.4e-4, sd=0.5, analyteConcentrations=c(1E-9, 2E-9, 1E-8, 2E-8, 1E-7), 
                        Rmax=250, associationLength=500, dissociationLength=500)
csm<-ConformationalSelectionModel(kon=1.72e5, koff=12.5e-2,kr=1.34e-1,kf=4e5, analyteConcentrations=c(1E-9, 2E-9, 1E-8, 2E-8, 1E-7), 
                                 Rmax=250, associationLength=500, dissociationLength=500)

lm_sim<-Simulate(lm, sampleFreq=0.5, sd=0.5)
ifm_sim<-Simulate(ifm, sampleFreq=0.5, sd=0.5)
csm_sim<-Simulate(csm, sampleFreq=0.5, sd=0.5)

dt_lm_sim<-lm_sim[[1]]

dt_ifm_sim1<-ifm_sim[[1]]
dt_ifm_sim2<-ifm_sim[[2]]
dt_ifm_sim<-GetObservedRUs(dt_ifm_sim1, dt_ifm_sim2)


dt_csm_sim1<-csm_sim[[1]]
dt_csm_sim2<-csm_sim[[2]]
dt_csm_sim<-GetObservedRUs(dt_csm_sim1, dt_csm_sim2)

plot(dt_lm_sim)
plot(dt_ifm_sim)
plot(dt_csm_sim)
