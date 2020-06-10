Date: June 8, 2020

*purpose: running pediatric HF analysis

*dependencies: flowCore, foreach, doParallel, iterators, plyr, randomForest, matrixStats, ROCR, FastKNN, miscTools, viridis, ggplot2, reshape2

* Change into your pediatric_HF directory as paths will be relative to this

*Note that I have this in script form as Run_PHF.R. Just modify the paths and the save directories and you can use that. 

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
```

# Step 2: Clustering Part=
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

Analysis of frequency differences. Again, change last argument to wherever you want to store this plot 

```R
source('VoPo_StandardAnalysis/vizFrequencyMap.R')
source('VoPo_StandardAnalysis/getFrequencyFeature.R')
FrF=getFrequencyFeature(Build,FNames)
uFrF=FrF[UnstimSamps,]
vizFrequencyMap(FrF,Layout,uResp,getwd())
```
