#####SPRTwoState package for analyzing and simulating SPR data
#  this R package is used to do the analysis under the general 
#  two model, especially for antigen and antibody interactions
#  
#  Version 0.1
#  For now, 1)it will only fit the unified two state model
#  2)It could also generalized the two state model for 
#    multiple conformations of antibodies
#  3) assume only antibodies (ligands) have multiple conformations, but
#     not antigens (analytes).
#
#### Developed By Feng @ BU 2016. All right reserved###


# @import MASS
#NULL

#roxygen2 comments below
#' @title S4 class definition of SensorgramData
#'
#' @description \code{SensorgramData} define the S4 class of SensorgramData
#'
#' @details defining the S4 class of the SensorgrameData.
#' This is the data structure to hold the sensor grame Data. 
#' It contains different slots for holding both association and 
#'    dissociation phase data.It also defines 
#'    the steady state window if there is one. Analyte concentrations 
#'    also included. Sometimes the data could 
#'    only have either association or dissociation phase alone. 
#'
#' @concept SPR
#'
#' @param associationData Dataframe holding the association phase data. 
#'         It has a fixed format/template.
#'         See \code{\link{ReadSensorgramData}} for details.
#' @param dissociationData Dataframe holding the dissociation phase data. 
#'          It has a fixed format/template.
#'          See \code{\link{ReadSensorgramData}} for details.
#' @param analyteConcentrations, numeric vector for analyteConcentrations 
#'          of unit molar concentration (M/L)
#' @param steadyStateStart numeric indicating the starting time of 
#'          steady state (optional).
#' @param steadyStateEnd numeric indicating the ending time of 
#'          steady state (optional).
#' @param offset numeric indicating the RU levels shift at the
#'          begining of dissociation phase. Not used for now.
#'          will implement in the future release. (optional)
#' @slot associationData data.frame
#' @slot dissociationData data.frame
#' @slot analyteConcentrations numeric vector
#' @slot steadyStateStart numeric
#' @slot steadyStateEnd numeric
#' @slot offset numeric
#' @seealso \code{\link{ReadSensorgramData}} \code{\link{GetObservedRUs}} \code{\link{SaveSPRData}}
#'
#' @examples
#'		dt<-new("SensorgramData",
#'		dissociationData=data.frame(time=1:5,RU=1:5, time2=1:5, RU=1:5))
#'		dt@associationData
#'
#'		dataFile<-list.files(system.file("extdata", package="SPRATS"),
#'		 pattern = "CBSInhibitor_raw.txt", full.names=TRUE) 
#'		rawData<-ReadSensorgramData(dataFile, skip=4, 
#'		header=T, sep="\t", associationPhaseEnd=80, dissociationPhaseEnd=150)
#'
#'		plot(rawData)
#'
#' @export
setClass( "SensorgramData",
		representation(associationData="data.frame",
					   dissociationData="data.frame",
					   analyteConcentrations="numeric",
					   steadyStateStart="numeric", #optional, but necessary for fitting two state
					   steadyStateEnd="numeric",#optional, but necessary for fitting two state
					   offset="numeric" ###optional, might be used in the future
		               ),
		prototype(associationData=data.frame(),
					   dissociationData=data.frame(),
					   analyteConcentrations=numeric(),
					   steadyStateStart=-1, 
					   steadyStateEnd=-1,
					   offset=0.0 
				)	   
		)

setValidity("SensorgramData", function(object){
				#first check for association data
				if(!is.null(object@associationData)&&length(object@associationData)!=0){
					#the columns must be pairs of Time and RUs
					if(dim(object@associationData)[2]%%2!=0)
					{
						cat("AssociationData is not in correct format!!\n")
						cat("dimension of association data:",dim(object@associationData)[2], "\n");
					}
					else if (!is.null(object@dissociationData)&&length(object@dissociationData)!=0){
					    if(dim(object@dissociationData)[2]%%2!=0)
							"DissociationData is not in correct format!!"
						else if(dim(object@associationData)[2]!=dim(object@dissociationData)[2]){
							"AssociationData and DissociationData are not consistent with each other!"
						}
						else {
							#cat("we are here\n")
							TRUE
						}
					}
					else 
						TRUE
				}
				#now check for dissociation data
				else if(!is.null(object@dissociationData)&&length(object@dissociationData)!=0){
					if(dim(object@dissociationData)[2]%%2!=0){
						"dissociationData is not in correct format!!"
					}
					else
						TRUE
				}
				#all the data is null or empty, go ahead
				
				else
				    TRUE
				#in this checking, we did not check for the number of concentration elements. leave it to the future
			}
		   )
		   
#' @title S3 method to read in the data
#'	
#' @description It read in data as a txt file according to
#'    exported SPR file from SensiQ
#' @details
#'  Data file has the following format\cr
#'	\tabular{llllllll}{
#'   Curve\tab	1\tab	2\tab		3\tab		4\cr	
#'   Analyte\tab	CBS\tab	CBS\tab		CBS\tab		CBS\cr	
#'  Ligand\tab	Ch 1\tab	Ch 1\tab		Ch 1\tab		Ch 1\cr	
#'  Conc\tab	0.00E+00\tab		5.00E-05\tab		0.00E+00\tab		0.00E+00\cr
#'	Time1\tab	Data1\tab	Time2\tab	Data2\tab	Time3\tab	Data3\tab	Time4\tab	Data4\cr
#'	-2.6\tab	0\tab	-2.6\tab	0.2839\tab	-2.6\tab	0\tab	-2.6\tab	0\cr
#'	-2.4\tab	0\tab	-2.4\tab	0.3632\tab	-2.4\tab	0\tab	-2.4\tab	0\cr
#'	-2.2\tab	0\tab	-2.2\tab	0.6978\tab	-2.2\tab	0\tab	-2.2\tab	0\cr
#'	}
#'	Each file contains the header about the information of the data
#'	Data colums for the raw sensorgram data with header for each column
#'	each "curve" has two column, one time and the other RUs. 
#'	The number of curves are variable. Please check the example
#'	data file in the example/ folder.

#'	Note: after reading the dissociation phase is rescaled to 
#'	start at time zero.!!! relative to the 
#'	dissociationPhaseStart/associationPhaseEnd
#'
#'	@param file the directory path to the input file with the
#'		correct format as mentioned above
#'	@param associationPhaseStart numeric indiciation the starting 
#'		time for associationPhase
#'	@param associationPhaseEnd numeric indiciation the starting 
#'		time for associationPhase. This is actually the starting
#'		time for dissociationPhase
#'	@param dissociationPhaseEnd numeric indiciation the starting 
#'		time for dissociationPhase
#' @param steadyStateStart numeric indicating the starting time of 
#'          steady state (optional).
#' @param steadyStateEnd numeric indicating the ending time of 
#'          steady state (optional).
#' @param analyteConcentrations numeric vector for analyte concentrations
#' @param skip numeric the options for control the text reading. indicating
#'		number of line to skip
#' @param header boolean header for the table reading (T/F)
#' @param sep character for the delimiter of the text
#'
#' @return a sensorgram data object
#'
#' @seealso  \code{\link{SensorgramData-class}} \code{\link{plot}} \code{\link{SaveSPRData}}
#' @export
ReadSensorgramData<-function(file, associationPhaseStart=0, associationPhaseEnd=-1,
								dissociationPhaseEnd=-1, steadyStateStart=-1, steadyStateEnd=-1,
								analyteConcentrations=numeric(0),
								skip=0, header=TRUE, sep="\t")
				{
					#first we need to check for validity of the input
					if(associationPhaseEnd<0 ||associationPhaseEnd< associationPhaseStart)
						cat("Warning: the associationPhase window is not correctly (?) set......\n")
					if(dissociationPhaseEnd<0||dissociationPhaseEnd<associationPhaseEnd)
						cat("Warning: the dissociationPhase window is not correctly (?) set......\n")
					cat("Start reading the data files.......\n")
					x_rawData<-read.table(file,skip=skip, header=header,sep=sep)
					cat("Done with reading the raw text file....\n")
					cat("Start parsing.....\n");
					x_sgdata<-new("SensorgramData")
					x_sgdata@analyteConcentrations<-analyteConcentrations
					
					#assuming formation is Time1 RU1 Time2 RU2....
					x_Ass<-data.frame();
					x_Dis<-data.frame();
					for(i in seq(1,dim(x_rawData)[2], by=2))
					{
						#cat("i:",i,"....\n")
						if(i+1>dim(x_rawData)[2])
							break;
							
						x_temp<-x_rawData[,c(i,i+1)]
						x_tempAss<-x_temp[x_temp[,1]>=associationPhaseStart&x_temp[,1]<=associationPhaseEnd,]
						x_tempDis<-x_temp[x_temp[,1]>=associationPhaseEnd&x_temp[,1]<=dissociationPhaseEnd,]
						x_tempDis[,1]<-x_tempDis[,1]-associationPhaseEnd
						if(i==1){
							x_Ass<-x_tempAss
							x_Dis<-x_tempDis
						}
						else
						{
							x_Ass<-cbind(x_Ass, x_tempAss)
							x_Dis<-cbind(x_Dis, x_tempDis)
						}
					}
					#now we have the data,
					cat("Setting sensorgram data......\n")
					x_sgdata@associationData<-x_Ass
					x_sgdata@dissociationData<-x_Dis
					x_sgdata@steadyStateStart<-steadyStateStart
					x_sgdata@steadyStateEnd<-steadyStateEnd
					cat("Done......\n")
					return(x_sgdata)
				}				

#' @title Generic plot function for SensorgramData
#' 
#' @description \code{plot} plot the sensorgram
#' @details Generic function to plot the sensorgram of SPR.
#'		Both the dissociation and association phase as well
#'		the steady state window.
#' 
#' @param SensorgramData the SensorgramData object to be plotted
#'
#' @seealso \code{\link{ReadSensorgramData}} \code{\link{SensorgramData-class}} \code{\link{SaveSPRData}}
#' @export
setMethod("plot", "SensorgramData",
		function(x, col=1, lty=1, lwd=1, main="",xlab="time (min)", ylab="RUs" )
		{
			#cat("times:",ceiling(length(x@associationData)/2/col))
			col<-rep(col,ceiling(length(x@associationData)/2/length(col)))
			lty<-rep(lty,ceiling(length(x@associationData)/2/length(lty)))
			lwd<-rep(lwd,ceiling(length(x@associationData)/2/length(lwd)))
			#cat("lty:",lty,"\n")
			max_x1<-0;
			max_x2<-0;
			max_y1<-0;
			max_y2<-0;
			meanRUs<-rep(0, length(x@associationData)/2)
			for(i in 1:(length(x@associationData)/2))
			{
				if(max_x1>min(x@associationData[,i*2-1]))
				{
					max_x1<-min(x@associationData[,i*2-1])
				}
				if(max_x2<max(x@associationData[,i*2-1]))
				{
					max_x2<-max(x@associationData[,i*2-1])
				}
				if(max_y1>min(x@associationData[,i*2]))
				{
					max_y1<-min(x@associationData[,i*2])
				}
				if(max_y2<max(x@associationData[,i*2]))
				{
					max_y2<-max(x@associationData[,i*2])
				}
				#meanRUs[i]<-
			}
			
			
			for(i in 1:(length(x@dissociationData)/2))
			{
				if(max_x1>min(max(x@associationData[,i*2-1])+x@dissociationData[,i*2-1]))
				{
					max_x1<-min(max(x@associationData[,i*2-1])+x@dissociationData[,i*2-1])
				}
				if(max_x2<max(max(x@associationData[,i*2-1])+x@dissociationData[,i*2-1]))
				{
					max_x2<-max(max(x@associationData[,i*2-1])+x@dissociationData[,i*2-1])
				}
				if(max_y1>min(x@dissociationData[,i*2]))
				{
					max_y1<-min(x@dissociationData[,i*2])
				}
				if(max_y2<max(x@dissociationData[,i*2]))
				{
					max_y2<-max(x@dissociationData[,i*2])
				}
				#meanRUs[i]<-
			}
			
			plot(c(max_x1, max_x2), c(max_y1, max_y2), type="n",main=main, xlab=xlab, ylab=ylab)
			
			for(i in 1:(length(x@associationData)/2))
			{
			#cat("plotting i:",i,"lty[i]:",lty[i],"\n")
				lines(x@associationData[,i*2-1], x@associationData[,i*2],col=col[i], lty=lty[i], lwd=lwd[i])
				lines(max(x@associationData[,i*2-1])+x@dissociationData[,i*2-1], x@dissociationData[,i*2],col=col[i],lty=lty[i],lwd=lwd[i])
			}
		}
	)
#'@title S3 function to combine/add two SensorgramData
#'
#'@description Combine two SensorgramData object to simulate
#'		the real readout on a SPR machine.  
#'
#'@details The idea is to combine the RUs from two different source
#'		e.g. from two complexes with different conformations, as
#'		an real readout on the SPR machine. It adds up the RUs at
#'		the same time points to get the sum. Do this through the 
#'		association and dissociation. Don't try to modify any other 
#'		slots for the object. 
#'
#'@param DATA1 SensorgramData 
#'@param DATA2 SensorgramData 

#'@return a SensorgramData object with the sum of the two input objects
#'		on both association and dissociation data. It takes other slot
#'		values from DATA1 assuming DATA2 has the identical values on
#'		these other slots. 
#'@seealso \code{\link{ReadSensorgramData}} \code{\link{SensorgramData-class}} \code{\link{plot}} \code{\link{SaveSPRData}}
#' @export
GetObservedRUs<-function(DATA1, DATA2)
	{
		if(class(DATA1)!="SensorgramData"||class(DATA2)!="SensorgramData")
		{
			cat("*******ERROR:SensorgramData are expected for both input")
			return(FALSE)
		}
		if(length(DATA1@associationData)!=length(DATA2@associationData))
		{
			cat("********ERROR:input data are not consistent with each other\n")
			return(FALSE);
		}
		if(length(DATA1@associationData[,1])!=length(DATA2@associationData[,1]))
		{
			cat("********ERROR:input data are not consistent with each other\n")
			return(FALSE);
		}
		#add all the AB together to get RU observed
		dt_RUs<-DATA1
		for(i in 1:length(DATA1@analyteConcentrations))
		{
			dt_RUs@associationData[,i*2]<-DATA1@associationData[,i*2]+DATA2@associationData[,i*2]
			dt_RUs@dissociationData[,i*2]<-DATA1@dissociationData[,i*2]+DATA2@dissociationData[,i*2]
		}
		dt_RUs
	}
#'@title S3 function to save SensorgramData data
#'
#'@description save SensorgramData object data to the disk.
#'		  
#'@details It save the object data to the disk in a text format,
#'		following the SensiQ format. See the detail \code{\link{ReadSensorgramData}}
#'		Note: it will issue a warning to complain about appending the table.
#'		There is nothing wrong with it. Simple ignore this warning.
#'
#'@param DATA SensorgramData to be saved 
#'@param file the directory path to the file
#'@param sep character to delimit the text fields 

#'
#'@seealso \code{\link{ReadSensorgramData}} \code{\link{SensorgramData-class}} \code{\link{plot}}	
#' @export
SaveSPRData<-function(DATA, file, sep="\t")
	{
		if(class(DATA)!="SensorgramData")
		{	
			cat("***ERROR***:sensorgramData is expected for the input\n")
			return(FALSE)
		}
		dtfm_RUs_a<-DATA@associationData
		dtfm_RUs_d<-DATA@dissociationData
		#check the dissociation data first, if the first row is zero for time, we need to combine this one 
		#with the last one of the association data to make it consistent
		if(dtfm_RUs_d[1,1]==0)
		{
			#get the mean of both values
			dtfm_RUs_a[length(dtfm_RUs_a[,1]),seq(2,length(dtfm_RUs_a),2)]<-
				(dtfm_RUs_a[length(dtfm_RUs_a[,1]),seq(2,length(dtfm_RUs_a),2)]+
					dtfm_RUs_d[1,seq(2,length(dtfm_RUs_d),2)])/2;
			#also remove the first one of the dissociation data
			dtfm_RUs_d<-dtfm_RUs_d[c(-1),]
		}
		#also need to add the time together for dissociation data
		for(i in 1:length(DATA@analyteConcentrations))
		{
			dtfm_RUs_d[,i*2-1]<-max(dtfm_RUs_a[,i*2-1])+dtfm_RUs_d[,i*2-1]
		}
		dtfm_RUs_a<-rbind(dtfm_RUs_a, dtfm_RUs_d)

		#now save it to the disk
		
		write(DATA@analyteConcentrations,file=file, sep=sep)
		write.table(dtfm_RUs_a, file=file, sep=sep,append=TRUE, row.names=FALSE, quote=FALSE)
		#TRUE
	}
	
				
				
				
#'@title S4 class of the SPR Langmuir model
#'
#'@description Langmuir model object
#'
#'@details check the manual for the detailed model specifications
#'		\url{http://}
#'
#'@slot kon numeric association rate constant
#'@slot koff numeric dissociation rate constant
#'@slot analyteConcentrations numeric vector for different
#'		analyte concentrations.
#'@slot s_d numeric standard deviation of the Gaussian noise. 
#'		by default s_d=0, which there is no noise.
#'@slot Rmax numeric the maximum RUs levels at the infinite
#'		high analyte concentration.
#'@slot associationLength numeric the time period to run the association
#'@slot dissociationLength numeric the time period to run the dissociation
#'@slot R0 numeric vector the starting RUs for the dissocaiton phase.
#'		optional
#'@slot offset the level shifted at the begining of the dissociation phase
#'		optional
#'@seealso \code{\link{ConformationalSelectionModel-class}}
#'		\code{\link{TwoStateModel-class}}
#'		\code{\link{InducedFitModel-class}}
#'
#' @export
setClass("LangmuirModel",
		representation(kon="numeric",
					   koff="numeric",
					   s_d="numeric",
					   analyteConcentrations="numeric", #vector
					   Rmax="numeric", 
					   associationLength="numeric", #option
					   dissociationLength="numeric", #option
					   R0="numeric",#vector optional, if not set the R0 will be getting values from the end of association phase
					   offset="numeric" ###optional, might be used in the future
		               ),
		prototype(kon=NULL,
					   koff=NULL,
					   s_d=0,
					   analyteConcentrations=numeric(), #vector
					   Rmax=NULL, 
					   associationLength=NULL,
					   dissociationLength=NULL,
					   R0=numeric(),#optional, if not set the R0 will be getting values from the end of association phase
					   offset=0 ###optional, might be used in the future
				)	   
		)
#'@title S4 generic help function to check model validity
#'
#'@description generic function to be called by constructor for 
#'		model validity
#'@details it actually is called by different model method 
#'		setValidty(objects) to check for the validity
#'
#'@param x the model object to be checked for
#'
#' @export 
setGeneric("CheckModelValidity", signature="object",
			function(object) standardGeneric("CheckModelValidity")
			)
#'@describeIn CheckModelValidity check the validity of LangmuirModel
#'
#'@export
setMethod("CheckModelValidity", c("object"="LangmuirModel"),
		function(object)
		{
			#first check for data integrity
				if(is.null(object@kon)){
					cat("kon is not set correctly\n")
					return(FALSE)
				}
				
				if(is.null(object@koff))
				{
					cat("koff is not set correctly\n")
					return(FALSE)
				}
		    #if(is.null(object@s_d))
		    #{
		    #  cat("Standard deviation is not set correctly\n")
		    #  return(FALSE)
		    #}
		    if(length(object@analyteConcentrations)<=0)
				{
					cat("analyteConcentrations are empty\n")
					return(FALSE)
				}
				if(sum(object@analyteConcentrations<0)>0)
				{
					cat("analyteConctrations are not set correctly\n")
					return(FALSE)
				}
				if(is.null(object@Rmax))
				{
					cat("Rmax is not set correctly\n")
					return(FALSE)
				}
				if(length(object@R0)>0&&sum(object@R0<0)>0)
				{
					cat("R0 are not set correctly\n")
					return(FALSE)
				}
				if(length(object@R0>0)&&length(object@analyteConcentrations)>0&&length(object@R0)!=length(object@analyteConcentrations))
				{
					cat("R0 and analyteConcentrations are not of equal length\n")
					return(FALSE)
				}
				#check the special case, where R0 has not been set, and we also did run association
				if(is.null(object@associationLength)&&is.null(object@dissociationLength)&&length(object@R0)==0)
				{
					cat("it has specified to run dissociation only, but R0 has not been set")
					return(FALSE)
				}
				#all others are optional.
				
				TRUE
		}
	)
setValidity("LangmuirModel",
			CheckModelValidity
			
		   )
#'@title Constructor for LangmuirModel
#'
#'@seealso \code{\link{LangmuirModel}}
#' @export
LangmuirModel<-function(kon,
						koff, s_d=0,
						analyteConcentrations,
						Rmax,#="numeric", 
					   associationLength,#="numeric", #option
					   dissociationLength,#="numeric", #option
					   R0=numeric(),#vector optional, if not set the R0 will be getting values from the end of association phase
					   offset=0 ###optional, might be used in the future
					   )
	{
						new("LangmuirModel", kon=kon, koff=koff, s_d=s_d, analyteConcentrations=analyteConcentrations,
							Rmax=Rmax,
							associationLength=associationLength, dissociationLength=dissociationLength,
							R0=R0, offset=offset)
	}

#induced fit
#'@title S4 class of the SPR InducedFit model
#'
#'@description InducedFit model object
#'
#'@details check the manual for the detailed model specifications
#'		\url{http://}
#'
#'@slot BaseModel a LangmuirModel object to hold the parameter
#'		see \code{\link{LangmuirModel}}
#'@slot kf numeric association rate constant
#'@slot kr numeric dissociation rate constant
#'@slot R02 numeric vector the starting RUs for the dissocaiton phase.
#'		optional
#'@seealso \code{\link{LangmuirModel}} 
#'		\code{\link{ConformationalSelectionModel-class}}
#'		\code{\link{TwoStateModel-class}}
#' @export
setClass("InducedFitModel",
		representation(BaseModel="LangmuirModel",
					   kf="numeric", #ka, forward, for conformation phase  
					   kr="numeric", #kd, reverse, for conformation phase
					   R02="numeric"
					   ),
		prototype(BaseModel=NULL,
					   kf=NULL,
					   kr=NULL,
					   R02=numeric()  #this R02 is used to for running dissociation, if we want the R0_AB_star, the stabler format of the complex
				)	   
		)

#'@describeIn CheckModelValidity check the validity of InducedFitModel
#'@export
setMethod("CheckModelValidity", c("object"="InducedFitModel"),
		function(object)
		{
			#first check for data integrity
				if(is.null(object@BaseModel))
				{
					return("BaseModel is not set correctly")
				}
				if(!CheckModelValidity(object@BaseModel))
					{
						return("BaseModel is not set correctly") 
					}
				
				if(is.null(object@kf)){
					return("kf is not set correctly")
					#return(FALSE)
				}
				
				if(is.null(object@kr))
				{
					return("koff is not set correctly")
					#return(FALSE)
				}
				if(length(object@R02)>0&&sum(object@R02<0)>0)
				{
					cat("R0 are not set correctly\n")
					return(FALSE)
				}
				if(length(object@R02>0)&&length(object@BaseModel@analyteConcentrations)>0&&length(object@R02)!=length(object@BaseModel@analyteConcentrations))
				{
					cat("R02 and analyteConcentrations are not of equal length\n")
					return(FALSE)
				}
				if(length(object@R02)!=length(object@BaseModel@R0))
				{
					cat("R02 and R01 are not of equal size")
					return(FALSE)
				}
				TRUE
		}
	)
setValidity("InducedFitModel",
			CheckModelValidity
		   )
#'@title Constructor for \code{InducedFitModel}
#'
#'@seealso \code{\link{InducedFitModel-class}}
#' @export		   
InducedFitModel<-function(kon,
				koff,kf, kr,s_d=0,
				analyteConcentrations,
				Rmax,#="numeric", 
			   associationLength,#="numeric", #option
			   dissociationLength,#="numeric", #option
			   R0=numeric(),#vector optional, if not set the R0 will be getting values from the end of association phase
			   R02=numeric(),
			   offset=0 ###optional, might be used in the future
			   )
	{
						lgm<-new("LangmuirModel", kon=kon,koff=koff, s_d=s_d, analyteConcentrations=analyteConcentrations,
							Rmax=Rmax,
							associationLength=associationLength, dissociationLength=dissociationLength,
							R0=R0, offset=offset)
						new("InducedFitModel", BaseModel=lgm, kf=kf, kr=kr, R02=R02)
	}
#conformational selection
#'@title S4 class of the SPR ConformationalSelection model
#'
#'@description ConformationalSelection model object
#'
#'@details check the manual for the detailed model specifications
#'		\url{http://}
#'
#'@slot BaseModel a LangmuirModel object to hold the parameter
#'		see \code{\link{LangmuirModel-class}}
#'@slot kf numeric association rate constant
#'@slot kr numeric dissociation rate constant
#'@slot R02 numeric vector the starting RUs for the dissocaiton phase.
#'		optional
#'@seealso \code{\link{LangmuirModel-class}} 
#'		\code{\link{InducedFitModel-class}}
#'		\code{\link{TwoStateModel-class}}
#' @export

setClass("ConformationalSelectionModel",
		representation(BaseModel="LangmuirModel",
					   kf="numeric", #ka, forward, for conformation phase  
					   kr="numeric" #kd, reverse, for conformation phase
					   ),
		prototype(BaseModel=NULL,
					   kf=NULL,
					   kr=NULL
				)	   
		)
#'@describeIn CheckModelValidity check the validity of ConformationalSelectionMdoel
#'@export
setMethod("CheckModelValidity", c("object"="ConformationalSelectionModel"),
		function(object)
		{
			#first check for data integrity
				if(is.null(object@BaseModel))
				{
					return("BaseModel is not set correctly")
				}
				if(!CheckModelValidity(object@BaseModel))
				{
					return("BaseModel is not set correctly")
				}
				if(is.null(object@kf)){
					return("kf is not set correctly")
					#return(FALSE)
				}
				
				if(is.null(object@kr))
				{
					return("koff is not set correctly")
					#return(FALSE)
				}
				
				TRUE
		}
	)
setValidity("ConformationalSelectionModel",
			CheckModelValidity
			) 
#'@title Constructor for \code{ConformationalSelectionModel}
#'
#'@seealso \code{\link{ConformatinalSelectionModel-class}}
#' @export		   
ConformationalSelectionModel<-function(kon,
				koff,kf, kr,s_d=0,
				analyteConcentrations,
				Rmax,#="numeric", 
			   associationLength,#="numeric", #option
			   dissociationLength,#="numeric", #option
			   R0=numeric(),#vector optional, if not set the R0 will be getting values from the end of association phase
			   
			   offset=0 ###optional, might be used in the future
			   )
	{
						lgm<-new("LangmuirModel", kon=kon,koff=koff, s_d=s_d, analyteConcentrations=analyteConcentrations,
							Rmax=Rmax,
							associationLength=associationLength, dissociationLength=dissociationLength,
							R0=R0, offset=offset)
						new("ConformationalSelectionModel", BaseModel=lgm, kf=kf, kr=kr)
	}		   
#two state model
#in this two state model case, we always using the parameters for
# length of period
# Rmax
# concentration
# R0
# offset
# from the inducedFitModel@LangmuirModel

#'@title S4 class of the SPR Two State model
#'
#'@description Two State model object
#'
#'@details check the manual for the detailed model specifications
#'		\url{http://}
#'		Note, 
#'		in this two state model case, we always using the parameters for
#'		length of period\cr
#'		 Rmax\cr
#'		 concentration\cr
#'		 R0\cr
#' 		offset\cr
#'		 from the inducedFitModel@LangmuirModel
#'
#'@slot BaseModel1 an InducedFitModel object to hold the parameter
#'		see \code{\link{InducedFitModel-class}}
#'@slot BaseModel2 a ConformationalSelectionModel object to hold the parameter
#'		see \code{\link{ConformationalSelectionModel-class}}

#'@seealso \code{\link{LangmuirModel-class}} 
#'		\code{\link{InducedFitModel-class}}
#'		\code{\link{ConformationalSelectionModel-class}}
#' @export

setClass("TwoStateModel",
		representation(BaseModel1="InducedFitModel",
					   BaseModel2="ConformationalSelectionModel"
					   ),
		prototype(BaseModel1=NULL,
				  BaseModel2=NULL
				)	   
		)
#'@describeIn CheckModelValidity check the validity of TwoStateModel
setMethod("CheckModelValidity", c("object"="TwoStateModel"),
		function(object)
		{
			#first check for data integrity
				if(is.null(object@BaseModel1))
				{
					return("BaseModel1 is not set correctly")
				}
				if(!CheckModelValidity(object@BaseModel1))
				{
					return("BaseModel1 is not set correctly")
				}
				if(is.null(object@BaseModel2))
				{
					return("BaseModel2 is not set correctly")
				}
				if(!CheckModelValidity(object@BaseModel2))
				{
					return("BaseModel2 is not set correctly")
				}				
				TRUE
		}
	)
setValidity("TwoStateModel",
			CheckModelValidity
		   ) 
#'@title Constructor for \code{TwoStateModel}
#'
#'@seealso \code{\link{TwoStateModel-class}}
#' @export		   
TwoStateModel<-function(kon, koff, kf, s_d=0, kr,#for induced fit model
				kon2,koff2, kf2,kr2, #for conformational selection model
				analyteConcentrations,
				Rmax,#="numeric", 
			   associationLength,#="numeric", #option
			   dissociationLength,#="numeric", #option
			   R0=numeric(),#vector optional, if not set the R0 will be getting values from the end of association phase
			   R02=numeric(),#vector optional, if not set the R02 will be getting values from the end of association phase
			   offset=0 ###optional, might be used in the future
			   )
	{
		lgm<-new("LangmuirModel", kon=kon,koff=koff, s_d=s_d, analyteConcentrations=analyteConcentrations,
			Rmax=Rmax,
			associationLength=associationLength, dissociationLength=dissociationLength,
			R0=R0, offset=offset)
		lgIF<-new("InducedFitModel", BaseModel=lgm, kf=kf, kr=kr, R02=R02)

		lgm2<-new("LangmuirModel", kon=kon2,koff=koff2, s_d=s_d, analyteConcentrations=analyteConcentrations,
			Rmax=Rmax,
			associationLength=associationLength, dissociationLength=dissociationLength,
			R0=R0, offset=offset)
		lgCS<-new("ConformationalSelectionModel", BaseModel=lgm2, kf=kf2, kr=kr2)
		
		new("TwoStateModel", BaseModel1=lgIF, BaseModel2=lgCS)
	}		 

#doing the simulation
#'@title S4 Generic method to Simulate SPR data
#'
#'@description Generic method to run the simulation based on the model
#'
#'@details each model has different biochemical dynamics to generate
#'		the data. The Langmuir, induced fit and conformational selection
#'		models simulate data analytically, but the two state model simulates
#'		through numerical integration. Please check the details here
#'		\url{http://}
#'
#'@param sampleFreq (numeric): the time frequency to collect the SPR data
#'
#'@return a list of \code{\link{SensorgramData-class}} data object 
#'		holding the SPR data for different component in the system
#'@seealso \code{\link{ConformatinalSelectionModel-class}}
#'			\code{\link{TwoStateModel-class}}
#'			\code{\link{InducedFitModel-class}}
#'			\code{\link{LangmuirModel-class}}
#'			\code{\link{SensorgramData-class}}
#' @export		   
setGeneric("Simulate", signature="x",
			function(x, sampleFreq=0.01, s_d=0, ...) standardGeneric("Simulate"))
#'@describeIn Simulate to simulate the SPR data based on a Langmuir Model
setMethod("Simulate", c("x"="LangmuirModel"),#, "sampleFreq"="numeric", "s_d"="numeric),
		function(x,sampleFreq=0.01,s_d=0)
		{##running the analytical solution of langmuirModel
			#check the data integrity 
			if(!CheckModelValidity(x))
			{
				return(NULL)
			}
			
			x_sgd<-new("SensorgramData")#for AB
			x_sgd_A<-new("SensorgramData")# for A
			#cat("Rmax:", x@Rmax)
			#x_time<-seq(0,x@associationLength,by=sampleFreq)
			#we also want to check for the validity
			if(x@associationLength>0) #we will do this association phase
			{
			
				x_Ass<-data.frame();
				x_Ass_A<-data.frame();
				x_time<-seq(0,x@associationLength,by=sampleFreq)
				noise<-rnorm(length(x_time),0,s_d)

				for(i in c(1:length(x@analyteConcentrations)))
				{
					temp<-x@Rmax*x@analyteConcentrations[i]/(x@koff/x@kon+x@analyteConcentrations[i])*(1-exp(-1*x_time*(x@kon*x@analyteConcentrations[i]+x@koff)))
					if(i==1)
					{
						x_Ass<-data.frame(Time=x_time, RU1=temp)
						x_Ass_A<-data.frame(Time=x_time, RU1=x@Rmax-temp)
					}
					else
					{
						x_Ass<-cbind(x_Ass,data.frame(Time=x_time, RU1=temp))
						x_Ass_A<-cbind(x_Ass_A,data.frame(Time=x_time, RU1=x@Rmax-temp))
					}
				}
				x_sgd@associationData<-x_Ass+noise
				x_sgd_A@associationData<-x_Ass_A
				
			}#end of association
			#cat("Simulating dissociation phase......\n")
			#print(x_sgd@associationData)
			if(x@dissociationLength>0)#we will do this dissociation phage
			{
				x_Diss<-data.frame();
				x_Diss_A<-data.frame();
				x_time<-seq(0,x@dissociationLength,by=sampleFreq)
				noise<-rnorm(length(x_time),0,s_d)
				
				for(i in c(1:length(x@analyteConcentrations)))
				{
					x_R0<-0;
					
					if(length(x@R0)!=0)
					{
						x_R0<-x@R0[i]
					}
					else
					{
						x_R0<-x_sgd@associationData[length(x_sgd@associationData[,1]),i*2]
						#cat("R0:", x_R0,"\n");
					}
					#R02 for RAB_star
					x_R0<-x_R0+x@offset
					
					temp<-x_R0*exp(-1*x_time*x@koff)

					if(i==1)
					{
						x_Diss<-data.frame(Time=x_time, RU1=temp)
						x_Diss_A<-data.frame(Time=x_time, RU1=x@Rmax-temp)
					}
					else
					{
						x_Diss<-cbind(x_Diss,data.frame(Time=x_time, RU1=temp))
						x_Diss_A<-cbind(x_Diss_A,data.frame(Time=x_time, RU1=x@Rmax-temp))
					}
				}
				
				x_sgd@dissociationData<-x_Diss+noise
				x_sgd_A@dissociationData<-x_Diss_A
			}
			x_sgd@analyteConcentrations<-x@analyteConcentrations
			x_sgd_A@analyteConcentrations<-x@analyteConcentrations

			return(list(R_AB=x_sgd,R_A=x_sgd_A))
		}#end of function
	)
#'@describeIn Simulate to simulate the SPR data based on a InducedFit Model
setMethod("Simulate", c("x"="InducedFitModel"),
		function(x,sampleFreq=0.01,s_d=0)
		{##running the analytical solution of langmuirModel
			#check the data integrity 
			if(!CheckModelValidity(x))
			{
				return(NULL)
			}
			
			#prepare for the calculations
			ka1<-x@BaseModel@kon
			kd1<-x@BaseModel@koff
			ka2<-x@kf
			kd2<-x@kr
			s_d<-x@BaseModel@s_d
			Rmax<-x@BaseModel@Rmax
			x_sgd<-new("SensorgramData") #unstable AB
			x_sgd_AB_star<-new("SensorgramData") #stable AB
			x_sgd_A<-new("SensorgramData") # A alone
			
			#x_time<-seq(0,x@associationLength,by=sampleFreq)
			#we also want to check for the validity
			if(x@BaseModel@associationLength>0) #we will do this association phase
			{
				x_Ass<-data.frame();#AB
				x_Ass_AB_star<-data.frame();
				x_Ass_A<-data.frame()
				x_time<-seq(0,x@BaseModel@associationLength,by=sampleFreq)
				noise<-rnorm(length(x_time),0,s_d)
				
				for(i in c(1:length(x@BaseModel@analyteConcentrations)))
				{
					conc<-x@BaseModel@analyteConcentrations[i]
					r1<-0.5*((ka1*conc+kd1+ka2+kd2)+sqrt((ka1*conc+kd1-(ka2+kd2))*(ka1*conc+kd1-(ka2+kd2))+4*kd1*ka2))
					r2<-0.5*((ka1*conc+kd1+ka2+kd2)-sqrt((ka1*conc+kd1-(ka2+kd2))*(ka1*conc+kd1-(ka2+kd2))+4*kd1*ka2))
					p<-ka1*conc*kd2+kd2*kd1+ka1*conc*ka2;
					m<-r1-r2;
					A<-ka1*conc*Rmax*((ka1*conc+kd1-r2)*(ka2+kd2)-ka2*kd1)/(p*m);
					B<-ka1*conc*Rmax*(kd1*ka2-(ka1*conc+kd1-r1)*(ka2+kd2))/(p*m);
  
					E_bar<-Rmax*kd1*kd2/p
					ES_star_bar<-Rmax*ka1*conc*ka2/p

					E<-E_bar+A*exp(-1*r1*x_time)+B*exp(-1*r2*x_time);
					ES_star<-ES_star_bar+1/kd1*(r1-(ka1*conc+kd1))*A*exp(-1*r1*x_time)+1/kd1*(r2-(ka1*conc+kd1))*B*exp(-1*r2*x_time);
					ES<-Rmax-E-ES_star

			
					if(i==1)
					{
						x_Ass<-data.frame(Time=x_time, RU1=ES)
						x_Ass_AB_star<-data.frame(Time=x_time, RU1=ES_star)
						x_Ass_A<-data.frame(Time=x_time, RU1=E)
					}
					else
					{
						x_Ass<-cbind(x_Ass,data.frame(Time=x_time, RU1=ES))
						x_Ass_AB_star<-cbind(x_Ass_AB_star,data.frame(Time=x_time, RU1=ES_star))
						x_Ass_A<-cbind(x_Ass_A,data.frame(Time=x_time, RU1=E))
					}
				}
				x_sgd@associationData<-x_Ass
				x_sgd_AB_star@associationData<-x_Ass_AB_star+noise
				x_sgd_A@associationData<-x_Ass_A
				
			}#end of association
			#cat("Simulating dissociation phase......\n")
			#print(x_sgd@associationData)
			if(x@BaseModel@dissociationLength>0)#we will do this dissociation phase
			{
				x_Diss<-data.frame();
				x_Diss_AB_star<-data.frame();
				x_Diss_A<-data.frame();
				x_time<-seq(0,x@BaseModel@dissociationLength,by=sampleFreq)
				noise<-rnorm(length(x_time),0,s_d)
				
				for(i in c(1:length(x@BaseModel@analyteConcentrations)))
				{
					x_R0<-0;
					X_R0_star<-0

					if(length(x@BaseModel@R0)!=0)
					{
						x_R0<-x@BaseModel@R0[i]
						x_R0_star<-x@R02
					}
					else
					{
						x_R0<-x_sgd@associationData[length(x_sgd@associationData[,1]),i*2]
						x_R02<-x_sgd_AB_star@associationData[length(x_sgd_AB_star@associationData[,1]),i*2]
						
						#cat("R0:", x_R0,"\n");
					}
					x_R0<-x_R0+x@BaseModel@offset
					x_R02<-x_R02+x@BaseModel@offset
					#=================
					conc0<-0;

					#calculate parameter again

					r1<--0.5*((ka1*conc0+kd1+ka2+kd2)+sqrt((ka1*conc0+kd1-(ka2+kd2))*(ka1*conc0+kd1-(ka2+kd2))+4*kd1*ka2))
					r2<--0.5*((ka1*conc0+kd1+ka2+kd2)-sqrt((ka1*conc0+kd1-(ka2+kd2))*(ka1*conc0+kd1-(ka2+kd2))+4*kd1*ka2))

					A<- -(kd1*x_R0+r2*(x_R0+x_R02))/(r2-r1)
					B<-  (kd1*x_R0+r1*(x_R0+x_R02))/(r2-r1)
					
					ES<-r1/kd1*A*exp(r1*x_time)+r2/kd1*B*exp(r2*x_time)
					ES_star<- -(r1+kd1)/kd1*A*exp(r1*x_time)-(r2+kd1)/kd1*B*exp(r2*x_time)
					E<-Rmax+A*exp(r1*x_time)+B*exp(r2*x_time)
					if(i==1)
					{
						x_Diss<-data.frame(Time=x_time, RU1=ES)
						x_Diss_AB_star<-data.frame(Time=x_time, RU1=ES_star)
						x_Diss_A<-data.frame(Time=x_time, RU1=E)
					}
					else
					{
						x_Diss<-cbind(x_Diss,data.frame(Time=x_time, RU1=ES))
						x_Diss_AB_star<-cbind(x_Diss_AB_star,data.frame(Time=x_time, RU1=ES_star))
						x_Diss_A<-cbind(x_Diss_A,data.frame(Time=x_time, RU1=E))
					}
				}
				
				x_sgd@dissociationData<-x_Diss
				x_sgd_AB_star@dissociationData<-x_Diss_AB_star+noise
				x_sgd_A@dissociationData<-x_Diss_A
			}
			x_sgd@analyteConcentrations<-x@BaseModel@analyteConcentrations
			x_sgd_AB_star@analyteConcentrations<-x@BaseModel@analyteConcentrations
			x_sgd_A@analyteConcentrations<-x@BaseModel@analyteConcentrations
			return(list(AB=x_sgd, AB_star=x_sgd_AB_star, A=x_sgd_A))
		}#end of function
	)
	
#the input is the model with parameters


#'@describeIn Simulate to simulate the SPR data based on a ConformationalSelection Model
#'
#'@details
#'the output is the sensorgram data for AB, A and A_star
#'   warning though, the association data for A and A_star 
#'   are good. The dissociation data for A is the sum of A and A_star
#'   since we can simply tell them apart with only R observed (R is the AB levels)
#'   So we simply add them together and save them in A@dissociationData and 
#'   save nothing in A_star (no dissociation data).
#'	Note: this part is a bit confusing. Need more work in the future.
setMethod("Simulate", c("x"="ConformationalSelectionModel"),
		function(x, sampleFreq=0.01, s_d=0)
		{
				#check the data integrity 
			if(!CheckModelValidity(x))
			{
				return(NULL)
			}
			
			#prepare for the calculations
			ka1<-x@BaseModel@kon
			kd1<-x@BaseModel@koff
			ka2<-x@kf
			kd2<-x@kr
			s_d<-x@BaseModel@s_d
			Rmax<-x@BaseModel@Rmax
			x_sgd<-new("SensorgramData") #stable AB
			x_sgd_A_star<-new("SensorgramData") #stable A
			x_sgd_A<-new("SensorgramData") # A alone
			
			#we also want to check for the validity
			if(x@BaseModel@associationLength>0) #we will do this association phase
			{
				x_Ass<-data.frame();#AB
				x_Ass_A_star<-data.frame();
				x_Ass_A<-data.frame();
				x_time<-seq(0,x@BaseModel@associationLength,by=sampleFreq)
				noise<-rnorm(length(x_time),0,s_d)
				
				for(i in c(1:length(x@BaseModel@analyteConcentrations)))
				{
					conc<-x@BaseModel@analyteConcentrations[i]
					r1<-0.5*((ka1*conc+kd1+ka2+kd2)+sqrt((ka1*conc+kd1-(ka2+kd2))*(ka1*conc+kd1-(ka2+kd2))+4*kd2*ka1*conc))
					r2<-0.5*((ka1*conc+kd1+ka2+kd2)-sqrt((ka1*conc+kd1-(ka2+kd2))*(ka1*conc+kd1-(ka2+kd2))+4*kd2*ka1*conc))
					p<-kd2*kd1+ka2*kd1+ka2*conc*ka1;
					m<-ka2+kd2;#r1-r2;
					h<-r1-r2;
					D1<--1*ka2*ka1*conc*kd2*Rmax*r2/(p*m*h);
					D2<-ka2*kd2*ka1*conc*Rmax*r1/(p*m*h);

					A_bar<-Rmax*kd1*kd2/p
					#A_star_bar<-Rmax*ka2*kd1/p
					AB_bar<-Rmax*ka1*conc*ka2/p;
					A<-A_bar+D1*exp(-1*r1*x_time)+D2*exp(-1*r2*x_time)
					AB<-AB_bar+1/kd2*(r1-(kd2+ka2))*D1*exp(-1*r1*x_time)+1/kd2*(r2-(kd2+ka2))*D2*exp(-1*r2*x_time);
					
					#cat("length of A:", length(A),"\n")
					#cat("length of AB:", length(AB),"\n")
					#cat("length of Rmax:", length(Rmax), "\n")
					
					A_star<-Rmax-A-AB;
					
					if(i==1)
					{
						x_Ass<-data.frame(Time=x_time, RU1=AB)
						x_Ass_A_star<-data.frame(Time=x_time, RU1=A_star)
						x_Ass_A<-data.frame(Time=x_time, RU1=A)
					}
					else
					{
						x_Ass<-cbind(x_Ass,data.frame(Time=x_time, RU1=AB))
						x_Ass_A_star<-cbind(x_Ass_A_star,data.frame(Time=x_time, RU1=A_star))
						x_Ass_A<-cbind(x_Ass_A,data.frame(Time=x_time, RU1=A))
					}
				}
				x_sgd@associationData<-x_Ass+noise
				x_sgd_A_star@associationData<-x_Ass_A_star
				x_sgd_A@associationData<-x_Ass_A
				
			}#end of association
			#cat("Simulating dissociation phase......\n")
			#print(x_sgd@associationData)
			if(x@BaseModel@dissociationLength>0)#we will do this dissociation phase
			{
				x_Diss<-data.frame();
				#x_Diss_A_star<-data.frame();
				x_Diss_A<-data.frame();
				x_time<-seq(0,x@BaseModel@dissociationLength,by=sampleFreq)
				noise<-rnorm(length(x_time),0,s_d)
				
				for(i in c(1:length(x@BaseModel@analyteConcentrations)))
				{
					x_R0<-0;
					X_R0_star<-0

					if(length(x@BaseModel@R0)!=0)
					{
						x_R0<-x@BaseModel@R0[i]
						#x_R0_star<-x@R02
					}
					else
					{
						x_R0<-x_sgd@associationData[length(x_sgd@associationData[,1]),i*2]
						#x_R02<-x_sgd_AB_star@associationData[length(x_sgd_AB_star@associationData[,1]),i*2]
						
						#cat("R0:", x_R0,"\n");
					}
					x_R0<-x_R0+x@BaseModel@offset
					#x_R02<-x_R02+x@BaseModel@offset
					#=================
					#conc0<-0;

					#calculate parameter again
					ES<-x_R0*exp(-kd1*x_time)
					E_total<-Rmax-ES
					if(i==1)
					{
						x_Diss<-data.frame(Time=x_time, RU1=ES)
						x_Diss_A<-data.frame(Time=x_time, RU1=E_total)
						#x_Diss_A_star<-data.frame(Time=x_time, RU1=E)
					}
					else
					{
						x_Diss<-cbind(x_Diss,data.frame(Time=x_time, RU1=ES))
						#x_Diss_A_star<-cbind(x_Diss,data.frame(Time=x_time, RU1=ES_star))
						x_Diss_A<-cbind(x_Diss_A,data.frame(Time=x_time, RU1=E_total))
					}
				}
				
				x_sgd@dissociationData<-x_Diss+noise
				x_sgd_A_star@dissociationData<-data.frame();
				x_sgd_A@dissociationData<-x_Diss_A
			}
			x_sgd@analyteConcentrations<-x@BaseModel@analyteConcentrations
			x_sgd_A_star@analyteConcentrations<-x@BaseModel@analyteConcentrations
			x_sgd_A@analyteConcentrations<-x@BaseModel@analyteConcentrations
			return(list(AB=x_sgd,  A=x_sgd_A, A_star=x_sgd_A_star))
		}
	)
#'@describeIn Simulate to simulate the SPR data based on a TwoState Model
setMethod("Simulate", c("x"="TwoStateModel"),
		function(x, sampleFreq=0.01, timeStep=0.01, s_d=0)
		{
			#check the data integrity 
			if(!CheckModelValidity(x))
			{
				return(NULL)
			}
			
			#prepare for the calculations
			  #IF section
			kon_i1<-x@BaseModel1@BaseModel@kon#1E4#1E5;
			kd_i1<- x@BaseModel1@BaseModel@koff#0.005#0.005;
			ka_i2<-x@BaseModel1@kf#0.001#0.003
			kd_i2<-x@BaseModel1@kr#0.0002;
			  #CS section
			kon_c1<-x@BaseModel2@BaseModel@kon#1E5#1E4
			kd_c1<-x@BaseModel2@BaseModel@koff#0.001;
			ka_c2<-x@BaseModel2@kf#0.0002;#0.002
			kd_c2<-x@BaseModel2@kr#0.004;
			
			Rmax<-x@BaseModel1@BaseModel@Rmax
			
			if(sampleFreq<timeStep)
			{
				sampleFreq<-timeStep;
			}
			
			x_sgd<-new("SensorgramData") #unstable AB
			x_sgd_AB_star<-new("SensorgramData")#stable AB
			x_sgd_A_star<-new("SensorgramData") #stable A
			x_sgd_A<-new("SensorgramData") # A alone
			#cat("kon_i1:",kon_i1,";kd_i1:", kd_i1,";ka_i2:", ka_i2,";kd_i2:", kd_i2,"\n");
			#cat("kon_c1:",kon_c1,";kd_c1:", kd_c1,";ka_c2:", ka_c2,";kd_c2:", kd_c2,"\n");
			#cat("Rmax:",Rmax,"\n");
			#we also want to check for the validity
			if(x@BaseModel1@BaseModel@associationLength>0) #we will do this association phase
			{
				cat("simulating the association phase data......\n");
				x_Ass<-data.frame();#AB
				x_Ass_AB_star<-data.frame();
				x_Ass_A_star<-data.frame();
				x_Ass_A<-data.frame();
				x_time<-seq(0,x@BaseModel1@BaseModel@associationLength,by=timeStep)
				noise<-rnorm(length(x_time),0,s_d)
			
				#cat("simulating assocation phase.......\n")
				for(i in c(1:length(x@BaseModel1@BaseModel@analyteConcentrations)))
				{
					conc<-x@BaseModel1@BaseModel@analyteConcentrations[i]

					#cat("i:", i,"\n")
					#cat("length of AB:", length(AB),"\n")
					#cat("length of Rmax:", length(Rmax), "\n")
					
					R_B_star<-x_time;
					R_B<-x_time;
					R_AB<-x_time;
					R_AB_star<-x_time;
					
					#Ratio_B<-kd_c2/ka_c2;
					R_B_star[1]<-Rmax/(kd_c2/ka_c2+1);
					R_B[1]<-Rmax - R_B_star[1];
					R_AB[1]<-0
					R_AB_star[1]<-0
					#R_B_prime[1]<-R_B[1];

				
					for( j in c(1:length(R_AB)))
					{
						R_AB_star_i<-(Rmax-R_B_star[j]-R_AB[j]-R_B[j])
						deltaR_AB<-kon_i1*conc*R_B[j]-kd_i1*R_AB[j]-ka_i2*R_AB[j]+kd_i2*R_AB_star_i;
						
						deltaR_B<- -1*kon_i1*conc*R_B[j]+kd_i1*R_AB[j]-ka_c2*R_B[j]+kd_c2*R_B_star[j];
						deltaR_B_star<- ka_c2*R_B[j]-kd_c2*R_B_star[j]-kon_c1*conc*R_B_star[j]+kd_c1*R_AB_star_i;
						
						if(j<=(length(R_AB)-1))
						{
							R_AB[j+1]<-deltaR_AB*timeStep+R_AB[j];
							R_B_star[j+1]<-deltaR_B_star*timeStep+R_B_star[j];
							R_B[j+1]<-deltaR_B*timeStep+R_B[j]
							R_AB_star[j+1]<-(Rmax-R_B_star[j+1]-R_AB[j+1]-R_B[j+1]);
						}
					}
#					cat("done.......\n");
#cat("length x_time:",length(x_time),"\n");					
#cat("length R_AB:",length(R_AB),"\n");
#cat("length R_aB_star:",length(R_AB_star),"\n");
#cat("length R_B:",length(R_B),"\n");
#cat("length R_B_star:",length(R_B_star),"\n");
					if(i==1)
					{
						#cat("1\n")
						x_Ass<-data.frame(Time=x_time, RU1=R_AB)
						x_Ass_AB_star<-data.frame(Time=x_time, RU1=R_AB_star)
						x_Ass_A_star<-data.frame(Time=x_time, RU1=R_B_star)
						x_Ass_A<-data.frame(Time=x_time, RU1=R_B)
						#cat("2\n")
					}
					else
					{
						#cat("4\n")
						x_Ass<-cbind(x_Ass,data.frame(Time=x_time, RU1=R_AB))
						x_Ass_AB_star<-cbind(x_Ass_AB_star,data.frame(Time=x_time, RU1=R_AB_star))
						x_Ass_A_star<-cbind(x_Ass_A_star,data.frame(Time=x_time, RU1=R_B_star))
						x_Ass_A<-cbind(x_Ass_A,data.frame(Time=x_time, RU1=R_B))
					}
					
				}

				outDataRows<-seq(1,length(x_Ass[,1]),by=sampleFreq/timeStep);
				
				x_sgd@associationData<-x_Ass[outDataRows,]
				x_sgd_AB_star@associationData<-x_Ass_AB_star[outDataRows,]
				x_sgd_A_star@associationData<-x_Ass_A_star[outDataRows,]
				x_sgd_A@associationData<-x_Ass_A[outDataRows,]
				
			}#end of association
			cat("Simulating dissociation phase......\n")
			#print(x_sgd@associationData)
			if(x@BaseModel1@BaseModel@dissociationLength>0)#we will do this dissociation phage
			{
				x_Diss<-data.frame();
				x_Diss_AB_star<-data.frame();
				x_Diss_A<-data.frame();
				x_Diss_A_star<-data.frame();
				x_time<-seq(0,x@BaseModel1@BaseModel@dissociationLength,by=timeStep)
				noise<-rnorm(length(x_time),0,s_d)
				
								
				for(i in c(1:length(x@BaseModel1@BaseModel@analyteConcentrations)))
				{

				  ru_detach_AB<-x_time;
					ru_detach_AB_star<-x_time
					ru_detach_AB[1]<-x_sgd@associationData[length(x_sgd@associationData[,1]),i*2];
					ru_detach_AB_star[1]<-x_sgd_AB_star@associationData[length(x_sgd_AB_star@associationData[,1]),i*2];

					ru_detach_B<-x_time
					ru_detach_B_star<-x_time
					ru_detach_B[1]<-x_sgd_A@associationData[length(x_sgd_A@associationData[,1]),i*2];
					ru_detach_B_star[1]<-x_sgd_A_star@associationData[length(x_sgd_A_star@associationData[,1]),i*2];

					conc0<-0;
					for( j in c(1:length(ru_detach_AB)))
					{
						ru_detach_AB_star_i<-(Rmax-ru_detach_B_star[j]-ru_detach_AB[j]-ru_detach_B[j])
						deltaR_AB<-kon_i1*conc0*ru_detach_B[j]-kd_i1*ru_detach_AB[j]-ka_i2*ru_detach_AB[j]+kd_i2*ru_detach_AB_star_i;
						
						deltaR_B<- -1*kon_i1*conc0*ru_detach_B[j]+kd_i1*ru_detach_AB[j]-ka_c2*ru_detach_B[j]+kd_c2*ru_detach_B_star[j];
						deltaR_B_star<- ka_c2*ru_detach_B[j]-kd_c2*ru_detach_B_star[j]-kon_c1*conc0*ru_detach_B_star[j]+kd_c1*ru_detach_AB_star_i;
						
						if(j<=(length(ru_detach_AB)-1))
						{
							ru_detach_AB[j+1]<-deltaR_AB*timeStep+ru_detach_AB[j];
							
							ru_detach_B_star[j+1]<-deltaR_B_star*timeStep+ru_detach_B_star[j];
							ru_detach_B[j+1]<-deltaR_B*timeStep+ru_detach_B[j]
							ru_detach_AB_star[j+1]<-(Rmax-ru_detach_B_star[j+1]-ru_detach_AB[j+1]-ru_detach_B[j+1]);
						}
		
					}

					if(i==1)
					{
						x_Diss<-data.frame(Time=x_time, RU1=ru_detach_AB)
						x_Diss_AB_star<-data.frame(Time=x_time, RU1=ru_detach_AB_star)
						x_Diss_A<-data.frame(Time=x_time, RU1=ru_detach_B)
						x_Diss_A_star<-data.frame(Time=x_time, RU1=ru_detach_B_star)
					}
					else
					{
						x_Diss<-cbind(x_Diss,data.frame(Time=x_time, RU1=ru_detach_AB))
						x_Diss_AB_star<-cbind(x_Diss_AB_star,data.frame(Time=x_time, RU1=ru_detach_AB_star))
						x_Diss_A<-cbind(x_Diss_A,data.frame(Time=x_time, RU1=ru_detach_B))
						x_Diss_A_star<-cbind(x_Diss_A_star,data.frame(Time=x_time, RU1=ru_detach_B_star))
					}
				}
				outDataRows<-seq(1,length(x_Diss[,1]),by=sampleFreq/timeStep);
				x_sgd@dissociationData<-x_Diss[outDataRows,]
				x_sgd_AB_star@dissociationData<-x_Diss_AB_star[outDataRows,]
				x_sgd_A_star@dissociationData<-x_Diss_A_star[outDataRows,];
				x_sgd_A@dissociationData<-x_Diss_A[outDataRows,];
			}
			x_sgd@analyteConcentrations<-x@BaseModel1@BaseModel@analyteConcentrations
			x_sgd_AB_star@analyteConcentrations<-x@BaseModel1@BaseModel@analyteConcentrations
			x_sgd_A_star@analyteConcentrations<-x@BaseModel1@BaseModel@analyteConcentrations
			x_sgd_A@analyteConcentrations<-x@BaseModel1@BaseModel@analyteConcentrations
			cat("done.......\n")
			
			return(list(AB=x_sgd, AB_star=x_sgd_AB_star, A=x_sgd_A, A_star=x_sgd_A_star))
		}

	####now we can do the fitting.
	#mode, indicating what model will be used 
	#       2, two state steady state fitting by linear regression
	#		1, two state  fitting by nonlinear
#'@title fit the sensorgrame data based on the two state model
#'
#'@description estimate the mean equilibrium dissociation constant based on the two state model
#'
#'@details It can also estimate the mean dissociation rate constants
#'		 at the steady state. It assumes the two state model, 
#'		but can be easily generilized to
#'		more than "two-conformation" cases. In the multi-conformation cases,
#'		the mean dissociation rate constant can only be approximated
#'		empirically using the data of a small window at the beginning of
#'		the dissociation phase. Check detail of the implementation here
#'		\url{http://}
#'		Note: this one is a S4 method. it is probably not necessary. We will
#' 		change it(??).
#'@param x SensorgramData containing the data to be fitted
#'@param mode Integer to selection which type of fitting to do
#'		1 to do nonlinear regression\cr
#'		2 to do linear regression.\cr
#'		nonlinear regression is better
#'@param type integer to indicate which type of analysis in dissociation
#' 		phase to be used. 
#'		1. for two state do fitting and then regression. 
#'		2. for multi-state approximation
#'@param steadyStateStart numeric the starting time for steady state
#'		optional and if provided, will overwrite the one in sensorgramdata
#'@param steadyStateEnd numeric the ending time for steady state
#'		optional and if provided, will overwrite the ones in sensorgramdata
#'@param windowSize numeric the time period used to do approximation to
#'		estimate the mean dissociaiton rate constant
#'@param init.association list of initial values of parameter to do the 
#'		non-linear regression. It has the following members
#'		list(Rmax=200, KD=1E-3)
#'@param init.dissociation list of initial values of parameter to do the 
#'		non-linear regression. It has the following members for type 1 analysis
#'		list(R1=100, r1=0.1) and for type 2 analysis 
#'		list(R1=20, R2=200,r1=-0.005,r2=-0.001)
#'@param control list of control elements for run non-linear regression
#'@param trace boolean to control whether show the trace of
#'		nonlinear regression.
#'@return a list of parameters estimated
#'@export
			#input:

#setGeneric("FitTwoStateSPR", signature="x",
#			function(x,mode=1,...) standardGeneric("FitTwoStateSPR"))

	#    mode: see above "generic" definition
	#    steadyStateStart and steadyStateEnd: will override the values specified in the dataSet
	#	start:used to set the starting points for the nls for estimating dissociation phase values
	#   type: used to indicating whether we are doing the bi-conformaiton settings
	#' @importFrom MASS rlm
	#for import rlm function from MASS package
####setMethod("FitTwoStateSPR", c("x"="SensorgramData"),#, "sampleFreq"="numeric"),
FitTwoStateSPR<-function(x, mode=1,type=1,steadyStateStart=-1,steadyStateEnd=-1, windowSize=10,
			init.association=NULL,init.dissociation=NULL, 
			control=list(maxiter = 500,tol = 1e-2, minFactor=1/1E10)
			,trace=F
			)
		{
		
			#we first get data ready
			#if(class(x)!="SensorgramData")
			#{
			#	cat("the input data is not the correct SensorgramData\n");
			#	return(FALSE);
			#}
			#get the steady state window ready
			ssStart<-x@steadyStateStart
			ssEnd<-x@steadyStateEnd
			#cat("ssStart:",ssStart, ";ssEnd:", ssEnd,"\n");
			
			if(steadyStateStart>0)
			{
				#cat("in hereJ:", steadyStateStart,"\n")
				#cat(
				ssStart<-steadyStateStart
			}
			if(steadyStateEnd>0)
				ssEnd<-steadyStateEnd
			#cat("ssStart:",ssStart, ";ssEnd:", ssEnd,"\n");
			
			if(ssStart<0||ssStart<min(x@associationData[,1]))
			{
				cat("******ERROR:steady state start point is not set correctly...\n");
				return(FALSE)
			}
			if(ssEnd<0||ssEnd>max(x@associationData[,1]))
			{
				cat("******ERROR:steady state end point is not set correctly...\n");
				return(FALSE)
			}
			
			if(mode!=1&&mode!=2)
			{
				cat("******WARNING:unknonw value for fitting mode, use \"mode=1\" instead...\n");
				mode<-1;
				#stop("\n");
			}
			
			#by now, we have set the window correctly
			#next, let's get mean of the values to do the fitting
			op<-par(ask=T, #need to confirm for showing picture
				#mfrow=c(2,2),#2x2 pictures on one plot
				pty="s" #square plotting
			)
			plot(x)
			
			meanSPR<-rep(0,length(x@associationData[1,])/2);
			for(i in 1:length(meanSPR))
			{
				meanSPR[i]<-mean(x@associationData[x@associationData[,i*2-1]>=ssStart&x@associationData[,i*2-1]<=ssEnd,i*2])
				
				lines(c(ssStart,ssEnd),c(meanSPR[i], meanSPR[i]), col=2,lwd=2)
			}
			#now we get the means, need to do linear 
			concs<-x@analyteConcentrations
			
			if(mode==2)
			{
				concs_inv<-1/concs
				#do linear regression to get the KA
				meanSPR_inv<-1/meanSPR
				lr<-lm(meanSPR_inv~concs_inv)
				Rmax<-1/lr[[1]][1]
				KD<-lr[[1]][2]*Rmax
				mainstr<-expression(paste("Linear Regression:", frac(1,"RUs"),
								frac("kD","Rmax"),"*",frac(1,group("[","B","]")),"+", frac(1,"Rmax"),sep=""))
				
				plot(concs_inv,meanSPR_inv,main=mainstr
					,#frac(1,group("[","B","]")+frac(1,Rmax)"), 
					xlab=expression(frac(1,group("[","B","]"))), 
					,ylab="1/RUs", type="p")
				e_RUs<-KD*(1/Rmax)*concs_inv+1/Rmax
				lines(concs_inv,e_RUs, lty=2, col=2, lwd=2)
				legend(min(concs_inv),max(meanSPR_inv),legend=c("fitted"),lty=c(2),col=c(2),lwd=c(2))
				
			}else #if(mode==1) the default one
			{
				if(is.null(init.association))
				{
					start.association=list(Rmax=2*max(meanSPR), KD=2*max(concs))
				}
				nlr<-nls(meanSPR~Rmax*concs/(KD+concs), start=init.association,
							control = control, 
							trace = TRUE)
				slr<-summary(nlr)
				Rmax<-slr[["parameters"]][1,1]
				KD<-slr[["parameters"]][2,1]
				mainstr<-expression(paste("Nonlinear Regression:RUs=" ,
								frac(paste("Rmax*",group("[","B","]"),sep=""),paste("kD*",group("[","B","]"),sep=""))
								,sep=""))
				
				plot(concs,meanSPR,main=mainstr
					,xlab=expression(group("[","B","]")), 
					,ylab="RUs", type="p")
				e_concs<-seq(min(concs),max(concs), by=(max(concs)-min(concs))/100)
				e_RUs<-Rmax*e_concs/(KD+e_concs)
				lines(e_concs,e_RUs, lty=2, col=2, lwd=2)
				legend(min(concs),max(meanSPR),legend=c("fitted"),lty=c(2),col=c(2),lwd=c(2))
				
			}
			#cat("printing it..........\n")
			#cat("KD:",KD)
			#print(lr)
			#summary(lr)
			##the intercept is 1/Rmax
			##the slope is 1/Rmax*1/Ka
			outList<-list(Rmax=Rmax, KD=KD, Concs=concs, MeanSPR=meanSPR)
			#cat(outList)
			#for dissociation phase, we have 
			
			cat("estimating data series #",i,"......\n")
			
			#maxRU<-max(x@dissociationData[,i*2])
			#RUs<-RUs+rnorm(length(RUs),0,0.001)
			Rrs<-c()
			R0s<-c()
			if(type==2)
			{
				
				for(i in 1:length(x@analyteConcentrations))
				{
					if(x@analyteConcentrations[i]<=0)
					{
						#we do the conc=0 series in this case
						next
					}
					#cat("running ", i,"\n");
					#doing the nls for each dissociation phase
					#startList<-
					if(is.null(init.dissociation))
					{
						init.dissociation<-list(R1=20, R2=200,r1=-0.005,r2=-0.001)
					}
					RUs<-x@dissociationData[,i*2]
					times<-x@dissociationData[,i*2-1]
					#rus<-RUs[,i]
					#tms<-times[,i]
					nlr<-nls(RUs~R1*exp(r1*times)+R2*exp(r2*times), start=init.dissociation,
							control = control, 
							trace = trace)
					slr<-summary(nlr)
					R1<-slr[["parameters"]][1,1]
					R2<-slr[["parameters"]][2,1]
					r1<-slr[["parameters"]][3,1]
					r2<-slr[["parameters"]][4,1]
					Rrs<-c(Rrs,(R1*r1+R2*r2))
					R0s<-c(R0s,(R1+R2))
				}
				
				#then doing regression to estimate the KA
				mainstr<-expression(paste("Linear Regression:",R[1],r[1],"+",R[2],r[2],"=-kd*R0" ,sep=""))
								
				plot(R0s, Rrs, xlab="R0", main=mainstr,
					ylab=expression(paste(R[1],r[1],"+",R[2],r[2],sep="")),
					type="p")
				slr_lm<-rlm(Rrs~R0s-1)
				outList$Rrs<-Rrs
				outList$R0s<-R0s
				outList$kd<-slr_lm[[1]]*-1
				e_RUs<-outList$kd*-1*R0s
				lines(R0s,e_RUs, lty=2, col=2, lwd=2)
				legend(min(R0s),max(e_RUs),legend=c("fitted"),lty=c(2),col=c(2),lwd=c(2))
				
			}
			else
			{
				#estimate the empirical rate of dissociation
				#cat("not implemented yet..........\n")
				for(i in 1:length(x@analyteConcentrations))
				{
					if(x@analyteConcentrations[i]<=0)
					{
						#we do the conc=0 series in this case
						next
					}
					times<-x@dissociationData[,i*2-1]
					RUs<-x@dissociationData[,i*2]
					tms<-times[times<windowSize]
					rus<-RUs[times<windowSize]
					
					#doing the nls for each dissociation phase
					#startList<-
					if(is.null(init.dissociation))
					{
						init.dissociation<-list(R1=20, r1=0.005)
					}
					nlr<-nls(rus~R1*exp(-r1*tms), start=init.dissociation,
							control = control, 
							trace = trace)
					slr<-summary(nlr)
					R1<-slr[["parameters"]][1,1]
					#R2<-slr[["parameters"]][2,1]
					r1<-slr[["parameters"]][2,1]
					#r2<-slr[["parameters"]][4,1]
					Rrs<-c(Rrs,(r1))
					R0s<-c(R0s,(R1))
				}
				kd<-mean(Rrs)
				outList$R0s<-R0s
				outList$kd<-kd
				
				mainstr<-paste("approximation of mean kd\n(window size=",windowSize," secs)",sep="")
								
				plot(c(1:length(R0s)),Rrs,  xlab="i", main=mainstr,
					ylab="kd",
					type="p")
				lines(c(1:length(R0s)),rep(kd, length(R0s)),lty=2,col=2, lwd=2)
				legend(1,max(R0s),legend=c("fitted"),lty=c(2),col=c(2),lwd=c(2))
			}
			par(op);
			outList
		}#end of the function
	
