Date: June 8, 2020

purpose: running pediatric HF analysis

dependencies: flowCore, foreach, doParallel, iterators, plyr, randomForest, matrixStats, ROCR, FastKNN, miscTools, viridis, ggplot2, reshape2

* Change into your pediatric_HF directory as paths will be relative to this

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
for(i in 1:length(sList)){
Stim=c(Stim,sList[[i]][5])
Val=grepl('H',sList[[i]][4])
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
''' 


# Step 2: Clustering Part=
* We will cluster all FCS files together and do splitting of samples later
* This should take 10-15 mins to run (depending on what happens on nalab2)

#In this example we will be doing our clustering on phenotypic markers
ToUse=PhenoInds

source('VoPo_StandardAnalysis/runRepMetaclust.R')

#for cluster-level visualization, I like to use 200 iterations. 50 clusters and 1000 numCPF are default parameters
Build=runRepMetaclust(200,50,FileNames,doCPF='specify',numCPF=1000,MN,ToUse,35)

##################################
#Basic plots
####################################
#Let's make the VoPo plots so that we can see what is important

#Phenotype related plots
source('VoPo_StandardAnalysis/vizClusterPhenotype.R')
#Layout=vizClusterPhenotype(Build,ToUse,'~/Clean_BClust/pediatric_HF/Phenotype')

#Let's do analysis of unstim samples first
UnstimSamps=which(Stim=='Unstim')

#let's extract function-based features
source('VoPo_StandardAnalysis/vizFunctionMaps.R')
source('VoPo_StandardAnalysis/getFunctionalFeature.R')

#Analysis of functional differences#

#extract functional features from VoPo object
fFeat=getFunctionalFeature(Build,FNames,FuncInds)
ufFeat=fFeat[UnstimSamps,]

#and class labels that correspond
uResp=Class[UnstimSamps]

#make maps
vizFunctionMaps(Layout,ufFeat,MN,FuncInds,uResp,'~/Clean_BClust/pediatric_HF/Func_unstim')

#Analysis of frequency differences#
source('VoPo_StandardAnalysis/vizFrequencyMap.R')
source('VoPo_StandardAnalysis/getFrequencyFeature.R')
FrF=getFrequencyFeature(Build,FNames)
uFrF=FrF[UnstimSamps,]

vizFrequencyMap(FrF,Layout,uResp,'~/Clean_BClust/pediatric_HF')
