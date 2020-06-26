This will show you how to make single-cell visualization from the VoPo clustering result.

# Step 1: We assume you have the VoPo Clustering result (`Build`) from the README file

We will do the rest of the stuff here with the VoPo output object (`Build`

Read in all of the stuff you need.

```R
library('Rtsne')
library('flowCore')
source('VoPo_StandardAnalysis/SampleCells.R')
source('VoPo_StandardAnalysis/vizAtlas_Function.R')
source('VoPo_StandardAnalysis/vizAtlas_Phenotype.R')
source('VoPo_StandardAnalysis/vizAtlas_Freq.R')
```

# Step 2: Dimensionality Reduction part

First we will get a sample of 30,000 cells across all sample files. Note I have already done this for you so I will leave these commented out. Inputs `FileNames`, `MN`, `ToUse` are what they were in the README. 

```R
#CellMat=SampleCells(FileNames,MN,ToUse,1000,30000)
#tRes=Rtsne(CellMat)$Y
```

Here you can read-in pre-processed versions of those data so layout will be fixed. 

```R
CellMat=readRDS('CellMat')
tRes=readRDS('tRes')
```

# Unstim Analysis

```R
UnstimSamps=which(Stim=='Unstim')
```

Analysis of functional differences

```R
#extract functional features from VoPo object

fFeat=getFunctionalFeature(Build,FNames,FuncInds)
ufFeat=fFeat[UnstimSamps,]

#and class labels that correspond

uResp=Class[UnstimSamps]
```
Make single-cell map for function specifying the directory `myDirec`

```R
vizAtlas_Function(CellMat,Build,ufFeat,uResp,ToUse,SampsToUse=NULL,35,tRes,myDirec,MN[FuncInds],50)
```
Make plots for phenotype and save to directory `PhenDirec`

```R
vizAtlas_Phenotype(CellMat,tRes,PhenDirec)
```
Get frequency based features

```R
FrF=getFrequencyFeature(Build,FNames)
uFrF=FrF[UnstimSamps,]
```
Make frequency maps and substitute in your `myDirec` for where you want these saved

```R
vizAtlas_Freq(CellMat,Build,uFrF,uResp,ToUse,SampsToUse=NULL,35,tRes,myDirec,MN[FuncInds],50)
```
# Stim Analysis

Create Stim list like we did last time to match samples

```R
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
```
Loop through all stims making a function and frequency map. Modify paths to directories. (See `dirNames`, which i update for each stim in the loop)

```R
Direcs=c('GMCSF','IFNa','IL','LPS')

for(i in 1:length(StimList)){
direc=Direcs[i]
uResp=Class[StimList[[i]][,1]]

#frequency features
GCSF_FrF=FrF[StimList[[i]][,2],]-FrF[StimList[[i]][,1],]

#frequency map
dirNames=paste('~/Clean_BClust/pediatric_HF/',direc,sep='')
vizAtlas_Freq(CellMat,Build,GCSF_FrF,uResp,ToUse,SampsToUse=NULL,35,tRes,dirNames,MN[FuncInds],50)

#function features
GCSF_FuF=fFeat[StimList[[i]][,2],]-fFeat[StimList[[i]][,1],]

#function map
vizAtlas_Function(CellMat,Build,GCSF_FuF,uResp,ToUse,SampsToUse=NULL,35,tRes,dirNames,MN[FuncInds],50)

```


