Date: June 8, 2020

* purpose: running pediatric HF analysis

* dependencies: flowCore, foreach, doParallel, iterators, plyr, randomForest, matrixStats, ROCR, FastKNN, miscTools, viridis, ggplot2, reshape2

* Change into your pediatric_HF directory as paths will be relative to this

* Note that I have this in script form as Run_PHF.R. Just modify the paths and the save directories and you can use that. 

```R
library('flowCore')
```

## Step 1:Pre-processing: finding data, separating functional and phenotypic markers, etc
We are going to do our analysis on mono-nuclear cells, which has granulocytes gated out

let's add a path to our files: we will input this to VoPo when we do clustering

```R
FileNames=list.files('~/stanleyn@stanford.edu/pediatric_HF/MNC','.fcs',full.names=TRUE)
```

I also create a short-hand name for book-keeping purposes (feature matrix rownames)

```R
FNames=list.files('~/stanleyn@stanford.edu/pediatric_HF/MNC','.fcs',full.names=FALSE)
```

We are going to exclude the internal controls from our analysis

```R
ICInds=which(grepl('IC',FNames))
FileNames=FileNames[-ICInds]
FNames=FNames[-ICInds]
```

Now we will get the metadata corresponding to each file from the file name. This will allow us to record the Stim, and Class for each FCS file. 

```R
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
```
Get your marker names, put them in a comprehensible human-understandable format

```R
frame=read.FCS(FileNames[1]) 
MN=pData(parameters(frame))[,2] 
```

let's use our annotations that we specified for whether each marker is functional or phenotypic

```R
markAnn=read.csv('MN_annotate.csv',header=TRUE,stringsAsFactors=FALSE)
PhenoInds=which(markAnn[,2]==1)
FuncInds=which(markAnn[,2]==2)
```

you can check that these indeed map to the right columns
```R
print(MN[PhenoInds])
print(MN[FuncInds])
```

# Step 2: Clustering Part

* We will cluster all FCS files together and do splitting of samples later
* This should take 10-15 mins to run (depending on what happens on nalab2)

In this example we will be doing our clustering on phenotypic markers. We specify that here and load the function

```R
ToUse=PhenoInds
source('VoPo_StandardAnalysis/runRepMetaclust.R')
```
for cluster-level visualization, I like to use 200 iterations. 50 clusters and 1000 numCPF are default parameters

```R
Build=runRepMetaclust(200,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)
```

# Step 3: Basic plots

Let's make the VoPo plots so that we can see what is important

First, Phenotype related plots. We are coloring each cluster by its phenotypic marker expresssion.  

```R
source('VoPo_StandardAnalysis/vizClusterPhenotype.R')
Layout=vizClusterPhenotype(Build,ToUse,'~/Clean_BClust/pediatric_HF/Phenotype')
```

Let's do analysis of unstim samples first. Grab only the files that correspond to the stims

```R
UnstimSamps=which(Stim=='Unstim')
```

Load code to extract features from our VoPo clustering result related to function

```R
source('VoPo_StandardAnalysis/vizFunctionMaps.R')
source('VoPo_StandardAnalysis/getFunctionalFeature.R')
```

Use those features to analyze functional differences

```R
fFeat=getFunctionalFeature(Build,FNames,FuncInds)
ufFeat=fFeat[UnstimSamps,]
```

Get class labels that correspond to unstim samples

```R
uResp=Class[UnstimSamps]
```

Make functional marker maps. Replace the last input with the path to where you want to store these plots (though you can write them to that example directory

```R
vizFunctionMaps(Layout,ufFeat,MN,FuncInds,uResp,'pediatric_HF/Func_unstim')
```

Analysis of frequency differences. Again, change last argument to wherever you want to store this plot. Currently I have it writing to your current working directory (which should be pediatric_HF folder)

```R
source('VoPo_StandardAnalysis/vizFrequencyMap.R')
source('VoPo_StandardAnalysis/getFrequencyFeature.R')
FrF=getFrequencyFeature(Build,FNames)
uFrF=FrF[UnstimSamps,]
vizFrequencyMap(FrF,Layout,uResp,getwd())
```

# Classification

You can use the feature matrices together for a classification task for this unstim comparison

let's make a joint set of features

```R
Joint=cbind(uFrF,ufFeat)
```

Your response variable is `uResp`. Use your favorite cross validation technique with features, `Joint` and response `uResp`.


## Saving The unstim data matrix and response vector to the `FeatMats` directory

Run the following code to save the data matrix and repsonse vector to the `FeatMats` directory. You will see below that we do the same and also save and rda file and response vector for each stimulation

```R
#save data matrix
save(Joint,file='FeatMats/unstim_dataMatrix.rda')
#save response vector (sample classes)
save(uResp,file='FeatMats/unstim_responseVector.rda')
```

# Analysis for each stim

* In general, for each stim, we subtract the unstim features for each patient from their stim features to build a new matrix

* The first step is to create a match matrix for indices so we are subtracting correctly.

* The output will be something called stim list, which for each stim the ith row will record index for unstimsample in column 1 and stim sample in column 2

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

* This means that in the ith element of your list, you will have unstim indices (column 1) and corresponding stim indices in column 2. 

* Always good to check if it's right. First stim is GMCSF

* Let's see if indices actually map to the correct thing

* It should pull the unstim sample and GCSF

```R
Test1=StimList[[1]][5,1]
Test2=StimList[[1]][5,2]
print(FNames[Test1])
```

* Now you can get the associated feature matrices for each stim by subtracting their features

* for example for frequency map for GCSF (which is the first element in our list)

* again, we are subtracting the stim value 

```R
GCSF_FrF=FrF[StimList[[1]][,2],]-FrF[StimList[[1]][,1],]
```

Then you can make the map for frequency

```R
vizFunctionMaps(Layout,GCSF_FrF,MN,FuncInds,uResp,'~/Clean_BClust/pediatric_HF/Func_unstim')
```

## Generating a feature matrix and visualization map for all of the Stims

* For each stimulation, you will be able to write its feature matrix and the corresponding response vector to the directory FeatMats.

* For example, let's say you want LPS feature matrix and response vector. Feature matrix will be in the dorectory `FeatMat` as `LPS_dataMatrix.rda`. The response vector will be `LPS_responseVector.rda`.

* Once you have run the below loop you will have generated this for each stimulation.

```R
Direcs=c('GMCSF','IFNa','IL','LPS')
for(i in 1:length(StimList)){

print('here')
direc=Direcs[i]
uResp=Class[StimList[[i]][,1]]

#frequency features
GCSF_FrF=FrF[StimList[[i]][,2],]-FrF[StimList[[i]][,1],]
#vizFrequencyMap(GCSF_FrF,Layout,uResp,direc)

GCSF_FuF=fFeat[StimList[[i]][,2],]-fFeat[StimList[[i]][,1],]
#vizFunctionMaps(Layout,GCSF_FuF,MN,FuncInds,uResp,direc)

#create the joint matrix and save
Joint_Stim=cbind(GCSF_FuF,GCSF_FrF)
FName_Matrix=paste('FeatMats/',names(StimList)[i],'_','dataMatrix','.rda',sep='')
save(Joint,file=FName_Matrix)

#save response vector (created above as uResp)
FName_Class=paste('FeatMats/',names(StimList)[i],'_','responseVector','.rda',sep='')
save(uResp,file=FName_Class)


}
```
