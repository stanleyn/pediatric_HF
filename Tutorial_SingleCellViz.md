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
Make plots for phenotype

```R
vizAtlas_Phenotype(CellMat,tRes,'~/shit/Phen_test')
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
