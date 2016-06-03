library(plyr)
library(ggplot2)
#library(knitr)
loadSTM<-function(x) {
 
  #x<-".//countfiles/met6.csv"
  counts<- read.csv(x, header=T,stringsAsFactors=FALSE) #open file
  counts<- counts[-1,-1] #remove gene barcodes and no match as they won't be used
  fcounts<- as.matrix (counts [grep("\\.1$",names(counts))] ) #forward
  rcounts<- as.matrix (counts [grep("\\.2$",names(counts))] ) #reverse
  tcounts<-fcounts+rcounts
  ncols=dim(tcounts)[2] #How many columns in CSV?
  if(grepl("_4r",x)){
  nmice=4
  }
   else if(grepl("_2r",x)){
  nmice=2
  }
  else{
  nmice=3
  }
  days=(ncols-1)/nmice #Remove 1 (the input), then divide by 3 (num replicates) to get number of days
  #MTG
  countarray<-array(NA,dim=c(nmice,days,length(tcounts))) #Raw counts for MOUSE, DAY, GENE
  print(paste("Loaded sample ",x,", which has reads for",days,"days"));
  
  for (d in 1:days ) { #Loop to make the days, 1st timepoint is day 4, etc.
    for(m in 1:nmice){ #Mice from 1 to 3
      countarray[m,d,]=tcounts[,((d-1)*nmice+m)]#Populate the 3D array
    }
    }
  
  genes=counts$gene #List of gene IDs
  input=tcounts[,ncol(tcounts)] #Raw counts for input column
  
  countarray<-replace(countarray,countarray==0,NA) #Set all 0s to NA so they are not used for further analysis
  
  inputsum=sum(input) #How many counts in total in the input
  inputratio=input/inputsum #What ratio of the input does each gene make up?
  inputincluded=inputratio>0.001 #List of booleans as to whether each gene is over the threshold
  
  filteredinput=input[inputincluded] #Filter input ratios
  filteredinputratio=filteredinput/sum(filteredinput)
  filteredgenes=genes[inputincluded] #... and gene names
  filteredarray=countarray[,,which(inputincluded==TRUE)] #...and raw counts
  countsums<- aaply(filteredarray,c(1,2), function(x) sum(x,na.rm=TRUE)) #Sum total for each timepoint
  ratioarray<-sweep(filteredarray, c(1,2), countsums, "/") #Calculate ratios for each gene for each timepoint

  #Now we need to model the variance for each sample & timepoint
  dims=dim(ratioarray) #standard array dimensions

pracvar<- array(NA,dims) #empty array to hold the variances

##typical values
volumecollected=20 #volume in ul collected from mouse
RBCsperul=1e7
parasitaemias=c(0.0007,0.007,0.05,0.15,0.19) #expected parasitaemia per day
efficiencyofdnaextraction=0.5#Yield from purification (somewhat fudged now)
volumeofresuspension=c(50,100,100,100,100) #volume DNA is resuspended in
volumeoftemplate=1 #volume of this resuspended DNA used in PCR

PCRvariances=loadOrGenerateCachedPCRVariances()

 
for (m in 1:dims[1]){
  for (d in 1:dims[2]){
    for (g in 1:dims[3]){
	###IN THE FUTURE THIS WHOLE SECTION WILL BE REPLACED WITH A LOOK UP TABLE FOR:
	# WHAT DAY OF INFECTION
	# HOW MANY READS
	# => WHAT VARIANCE?
    #Calculate estimated number of parasite genomes that make it to PCR template:
    ParasiteGenomesInSample=volumecollected*RBCsperul*parasitaemias[d]*efficiencyofdnaextraction*(volumeoftemplate/volumeofresuspension[d])
    probabilityThisParasite=ratioarray[m,d,g] #assume probability of this barcode is the measured ratio in the population
	
    ExpectedNumThisParasite=ParasiteGenomesInSample*probabilityThisParasite #np, binomial. How many of this barcode we expect make it to PCR tube
    VarianceThisParasite=ParasiteGenomesInSample*probabilityThisParasite*(1-probabilityThisParasite) #np(1-p). The variance in the number above
    VarianceThisRatio=VarianceThisParasite*(1/(ParasiteGenomesInSample**2)) # laws of variance, multiplying by a constant to calculate the variance of the barcode ratio.
   
    
    PCRnormvar=PCRvariances[which.min(abs(PCRvariances$initmolecules - ExpectedNumThisParasite)),]$normalisedvar
    
    
    SequencingN=(countsums[m,d]/2)
    SequencingP=ratioarray[m,d,g]
    SequencingVar=SequencingN*SequencingP*(1-SequencingP)
    SequencingMean=SequencingN*SequencingP
    SequencingNormalisedMean=SequencingMean/SequencingMean # 1
    SequencingNormalisedVar=SequencingVar/(SequencingMean**2)
    CombinedVarianceThisRatio=VarianceThisRatio*SequencingNormalisedVar + VarianceThisRatio*(SequencingNormalisedMean**2) + SequencingNormalisedVar*(probabilityThisParasite**2)
    FullCombinedVarianceThisRatio=CombinedVarianceThisRatio*PCRnormvar + CombinedVarianceThisRatio + PCRnormvar*(probabilityThisParasite**2)
 
    pracvar[m,d,g]=VarianceThisRatio
       }
     }
  
}


relfitness<- array(NA,c(dims[1],dims[2]-1,dims[3]))
relfitnessvar<- array(NA,dim(relfitness))
#Now calculate relative fitnesses (genes relative to each other at each timepoint) and their variances
for (m in 1:dims[1]){
  for (d in 1:(dims[2]-1)){
    for (g in 1:dims[3]){
      
     tomorrow=ratioarray[m,d+1,g]
     today=ratioarray[m,d,g]
     vartomorrow=pracvar[m,d+1,g]
     vartoday=pracvar[m,d,g]
     relfitness[m,d,g]=tomorrow/today
     relfitnessvar[m,d,g]=(tomorrow**2/today**2)*(vartomorrow/tomorrow**2+vartoday/today**2)
    }
  }
  
}

gaussianMeanAndVariance<- function(vals, variances){
  df<-data.frame(value=vals,variance=variances)
  df<-df[complete.cases(df),]
  
  vals=df$value
  variances=df$variance
  
  precs=1/variances
  
  mean=sum(vals*precs)/sum(precs)
	
	var1=(1/sum(precs))*(1/(length(vals)-1))*sum((vals-mean)**2/variances)
  var2=1/sum(precs)
 # var1=0
 if(is.na(var1)){var1<-0}
  var<-max(var1,var2)
  return(c(mean,var))
}


#Now calculate normalised fitnesses and their variances
controlgenes=c("PBANKA_103780","PBANKA_051500","p230p-tag","PBANKA_051490","PBANKA_140160","PBANKA_110420","PBANKA_103440")
calibratedfitnesses=c(1,1,1,1,0.55,0.60,0.65)

#controlgenes=c("PBANKA_103780","PBANKA_051500","p230p-tag","PBANKA_051490")
#calibratedfitnesses=c(1,1,1,1) #Uncomment this to only use normal growth controls


controlgeneids=  which(filteredgenes %in% controlgenes==TRUE) #This is a list of the locations of the controls in the geneset
other<-match(filteredgenes[controlgeneids],controlgenes) #This lists, in the same order, the locations in the controlgeneset





absfitness<- array(NA,dim(relfitness))
absfitnessvar<- array(NA,dim(relfitness))

controlcalibratedarray<-array(NA,c(dims[1],dims[2],length(controlgeneids)))
controlcalibratedvararray<-array(NA,c(dims[1],dims[2],length(controlgeneids)))
controlnames<-controlgenes[other]
for (m in 1:dims[1]){
  for (d in 1:(dims[2]-1)){
	timepointcalibratedfitnesses=relfitness[m,d,controlgeneids]/calibratedfitnesses[other]
	timepointcalibratedfitnessesvar=relfitnessvar[m,d,controlgeneids]/(calibratedfitnesses[other]**2)
    controlcalibratedarray[m,d,]<-timepointcalibratedfitnesses
	controlcalibratedvararray[m,d,]<-timepointcalibratedfitnessesvar
	
	
    gauss=gaussianMeanAndVariance(timepointcalibratedfitnesses,timepointcalibratedfitnessesvar)
    controlmean=gauss[1]
    controlvar=gauss[2]
    for (g in 1:dims[3]){
      absfitness[m,d,g]=relfitness[m,d,g]/controlmean
      absfitnessvar[m,d,g]=(relfitness[m,d,g]**2/controlmean**2)*((relfitnessvar[m,d,g]/relfitness[m,d,g]**2) + (controlvar/controlmean**2))
    }
  }
  
}
#Now calculate a single value per gene (assuming no time effect)
singleabsfitness <- rep(NA, dims[3])
singleabsfitnessvar <- rep(NA, dims[3])

d6toinput <- rep(NA, dims[3])

for (g in 1:dims[3]){
  fitnessmatrix<-absfitness[,,g]
  variancematrix<-absfitnessvar[,,g]
  fitnesses=c(fitnessmatrix)
  variances=c(variancematrix)
  gauss=gaussianMeanAndVariance(fitnesses,variances)
  singleabsfitness[g]=gauss[1]
  singleabsfitnessvar[g]=gauss[2]
  d6toinput[g]=mean(ratioarray[,3,g])/filteredinputratio[g]
  
}
bcdk=which(filteredgenes=="PBANKA_103780")
if("PBANKA_103780"%in% filteredgenes){
normd6toinputA=d6toinput/d6toinput[bcdk]
}
else{
normd6toinputA <- rep(NA, dims[3])
}
bcdk=which(filteredgenes=="PBANKA_110420")
if("PBANKA_110420"%in% filteredgenes){
normd6toinputB=d6toinput/d6toinput[bcdk]
}
else{
normd6toinputB <- rep(NA, dims[3])
}

bcdk=which(filteredgenes=="PBANKA_051500")
if("PBANKA_051500"%in% filteredgenes){
normd6toinputC=d6toinput/d6toinput[bcdk]
}
else{
normd6toinputC <- rep(NA, dims[3])
}
absfitnessdf=data.frame(gene=filteredgenes,fitness=singleabsfitness,variance=singleabsfitnessvar,d6toinput=d6toinput,normd6toinputA=normd6toinputA,normd6toinputB=normd6toinputB,normd6toinputC=normd6toinputC,day4abs=ratioarray[1,1,],day5absmax=pmax(ratioarray[1,2,],ratioarray[2,2,],ratioarray[3,2,]),day6abs=ratioarray[1,3,])
absfitnessdf$lower=absfitnessdf$fitness-sqrt(absfitnessdf$variance)*2
absfitnessdf$upper=absfitnessdf$fitness+sqrt(absfitnessdf$variance)*2


dayabsfitness<-array(NA,dim=c(dims[2]-1,dims[3]))
dayabsfitnessvar<-array(NA,dim=c(dims[2]-1,dims[3]))
for (d in 1:(dims[2]-1)){
  
  
  
  for (g in 1:dims[3]){
    gauss=gaussianMeanAndVariance(absfitness[,d,g],absfitnessvar[,d,g])
    dayabsfitness[d,g]=gauss[1]
    dayabsfitnessvar[d,g]=gauss[2]
     }
}
  
  merge<-absfitnessdf;
 # colnames(merge)[which(colnames(merge)=="fitness")]="Relative.Growth.Rate"
  merge$file=x 
  filename<-tail(unlist(strsplit(as.character(x), "/")),1)
   filename<-head(unlist(strsplit(as.character(filename), "\\.")),1)
  merge$experiment=filename 
  merge$variance=merge$variance*3.2 # hack to try to be better
 arrays<-list(genes=filteredgenes,input=filteredinput,counts=filteredarray,ratios=ratioarray,ratiosvar=pracvar,absfitness=absfitness,absfitnessvar=absfitnessvar,controlarray=controlcalibratedarray,controlvararray=controlcalibratedvararray,controlnames=controlnames)
  #Find integration effiencies
 return(return(list(table=merge,extra=arrays)))
}


loadOrGenerateCachedPCRVariances<-function(){
  if(file.exists("./otherdata/pcrvariances.csv")){
    df2<-read.csv("./otherdata/pcrvariances.csv")
    
    print("Using cached PCR variances")
    return(df2)
  }
  else{
    print( "Regenerating PCR variances")
rbinom_safe <- function(n,size,prob,max.size=2^31) {
  maxlen <- max(length(size),length(prob),n)
  prob <- rep(prob,length.out=maxlen)
  size <- rep(size,length.out=maxlen)
  res <- numeric(n)
  bigvals <- size>max.size
  if (nbig <- sum(bigvals>0)) {
    m <- (size*prob)[bigvals]
    sd <- sqrt(size*prob*(1-prob))[bigvals]
    res[bigvals] <- round(rnorm(nbig,mean=m,sd=sd))
  }
  if (nbig<n) {
    res[!bigvals] <- rbinom(n-nbig,size[!bigvals],prob[!bigvals])
  }
  return(res)
}

simulatePCR<-function(initialmolecules,cycles,probability) {
  if(cycles==0){return(initialmolecules)}
  else{
    newmolecules=initialmolecules+rbinom_safe(1,initialmolecules,probability);
    return (simulatePCR(round(newmolecules),cycles-1,probability))
  }
}

numcycles=45
efficiency=0.8
geometricstepsize=1.3
rangetoconsider=30
replicates=500
df=data.frame(initmolecules=round(geometricstepsize^(0:(rangetoconsider-1))),normalisedvar=NA)
for (i in 1:rangetoconsider){
  finalmolecules=replicate(replicates,simulatePCR(df[i,1],numcycles,efficiency))
  variance=sd(finalmolecules)**2
  normalisedvariance=variance/(mean(finalmolecules)**2)
  df[i,2]=normalisedvariance
  #print(paste(molecules,normalisedvariance))
}

df2=aggregate(normalisedvar ~ initmolecules, data = df, mean)
write.csv(df2,"pcrvariances.csv")
return(df2)
}
}