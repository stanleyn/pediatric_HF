#Date: June 8, 2020
#purpose: running pediatric HF analysis

library('flowCore')

######################################
#Step 1:Pre-processing: finding data, separating functional and phenotypic markers, etc
######################################

#We are going to do our analysis on mono-nuclear cells, which has granulocytes gated out

#let's add a path to our files: we will input this to VoPo when we do clustering
FileNames=list.files('~/stanleyn@stanford.edu/pediatric_HF/MNC','.fcs',full.names=TRUE)

#I also create a short-hand name for book-keeping purposes (feature matrix rownames)
FNames=list.files('~/stanleyn@stanford.edu/pediatric_HF/MNC','.fcs',full.names=FALSE)

#We are going to exclude the internal controls from our analysis
ICInds=which(grepl('IC',FNames))
FileNames=FileNames[-ICInds]
FNames=FNames[-ICInds]

#now we will get the metadata corresponding to each file from the file name
sList=strsplit(FNames,'_')
Stim=c()
Class=c()
PID=c()
for(i in 1:length(sList)){
Stim=c(Stim,sList[[i]][5])
Val=grepl('H',sList[[i]][4])
PID=c(PID,sList[[i]][4])
if(Val==TRUE){
Class=c(Class,0)
}
else{Class=c(Class,1)}
}

#get your marker names, put them in a comprehensible human-understandable format
frame=read.FCS(FileNames[1]) 
MN=pData(parameters(frame))[,2] 

#let's use our annotations that we specified for whether each marker is functional or phenotypic
markAnn=read.csv('MN_annotate.csv',header=TRUE,stringsAsFactors=FALSE)

PhenoInds=which(markAnn[,2]==1)
FuncInds=which(markAnn[,2]==2)

#you can check that these indeed map to the right columns
print(MN[PhenoInds])
print(MN[FuncInds])

########################################
#Step 2: Clustering Part
#########################################
# We will cluster all FCS files together and do splitting of samples later
# This should take 10-15 mins to run (depending on what happens on nalab2)

#In this example we will be doing our clustering on phenotypic markers
ToUse=PhenoInds

source('VoPo_StandardAnalysis/runRepMetaclust.R')

#for cluster-level visualization, I like to use 200 iterations. 50 clusters and 1000 numCPF are default parameters
#Build=runRepMetaclust(200,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)

##################################
#Basic plots
####################################
#Let's make the VoPo plots so that we can see what is important

#Phenotype related plots
#source('VoPo_StandardAnalysis/vizClusterPhenotype.R')
#Layout=vizClusterPhenotype(Build,ToUse,'~/Clean_BClust/pediatric_HF/Phenotype')

#Let's do analysis of unstim samples first
#UnstimSamps=which(Stim=='Unstim')

#let's extract function-based features
#source('VoPo_StandardAnalysis/vizFunctionMaps.R')
#source('VoPo_StandardAnalysis/getFunctionalFeature.R')

#Analysis of functional differences#

#extract functional features from VoPo object
#fFeat=getFunctionalFeature(Build,FNames,FuncInds)
#ufFeat=fFeat[UnstimSamps,]

#and class labels that correspond
#uResp=Class[UnstimSamps]

#make maps
#vizFunctionMaps(Layout,ufFeat,MN,FuncInds,uResp,'~/Clean_BClust/pediatric_HF/Func_unstim')

#Analysis of frequency differences#
#source('VoPo_StandardAnalysis/vizFrequencyMap.R')
#source('VoPo_StandardAnalysis/getFrequencyFeature.R')
#FrF=getFrequencyFeature(Build,FNames)
#uFrF=FrF[UnstimSamps,]

#vizFrequencyMap(FrF,Layout,uResp,'~/Clean_BClust/pediatric_HF')

#stop('') #just see if you can produce to here

#############################
#classification
##############################

#You can use the feature matrices together for a classification task for this unstim comparison

#let's make a joint set of features
#Joint=cbind(FrF,uFrF)

#Your response variable is uResp. Use your favorite cross validation technique with features, Joint and response uResp.

###########################
#stim analyses
##########################

#In general, for each stim, we subtract the unstim features for each patient from their stim features to build a new matrix

#The first step is to create a match matrix for indices so we are subtracting correctly.
#The output will be something called stim list, which for each stim the ith row will record index for unstimsample in column 1 and stim sample in column 2

UStim=unique(Stim)[1:4] #there are 4 stims and we exclude unstim
PID=paste('p',PID,sep='_')
UPat=unique(PID)
StimList=list()

for(u in 1:length(UStim)){
StimCand=grep(UStim[u],FNames)
UnstimInd=c()
MatchInds=c()
for(j in 1:length(UPat)){
c1=grep('Unstim',FNames)
c2=which(PID==UPat[j])          
intVal=intersect(c1,c2)
if(length(intVal)==0){
UnStimInd=c(UnstimInd,NA)
}
else{UnstimInd=c(UnstimInd,intVal)}

c11=grep(UStim[u],FNames)
intVal2=intersect(c11,c2)
if(length(intVal2)==0){
MatchInds=c(MatchInds,NA)
}
else{MatchInds=c(MatchInds,intVal2)}
}

newArray=cbind(UnstimInd,MatchInds)
StimList[[u]]=newArray
}

names(StimList)=UStim


#This means that in the ith element of your list, you will have unstim indices (column 1) and corresponding stim indices in column 2. 

#Always good to check if it's right. First stim is GMCSF
#Let's see if indices actually map to the correct thing
#It should pull the unstim sample and GCSF
Test1=StimList[[1]][5,1]
Test2=StimList[[1]][5,2]
print(FNames[Test1])
print(FNames[Test2])

#created directories to store these plots
Direcs=c('GMCSF','IFNa','IL','LPS')
print('test')
for(i in 1:length(StimList)){

print('here')
direc=Direcs[i]
uResp=Class[StimList[[i]][,1]]

#frequency features
GCSF_FrF=FrF[StimList[[i]][,2],]-FrF[StimList[[i]][,1],]
vizFrequencyMap(FrF,Layout,uResp,direc)

GCSF_FuF=fFeat[StimList[[i]][,2],]-fFeat[StimList[[i]][,1],]
vizFunctionMaps(Layout,GCSF_FuF,MN,FuncInds,uResp,direc)


}



