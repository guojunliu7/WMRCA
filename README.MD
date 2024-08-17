# WMRCA Tutorial

## Introduction

The R package `WMRCA` integrates the current common computational biology methods for cancer subtypes identification and provides a standardized framework for cancer subtype analysis based multi-omics data, such as gene expression, miRNA expression, DNA methylation and others. See `vignette('WMRCA')` for the comprehensive user guide.

## Installation

To install the WMRCA package in R, simply type

   ```
   install.packages('devtools')
   library(devtools)
   install_github('guojunliu7/WMRCA')
   ```
if the installation invokes i386 arch, try to use:
   ```
   install_github('guojunliu7/WMRCA', INSTALL_opts="--no-multiarch")
   ```
To make the WMRCA package fully functional (especially the Open Babel and rjava related functionalities), we recommend the users some tips: 

[>>Common installing tips](http://projects.scbdd.com/si/install_biomedr_tips.pdf)

Several dependencies of the WMRCA package may require some system-level libraries, check the corresponding manuals of these packages for detailed installation guides.
## Examples
Code snippet: example for mRNA data preprocessing
 ```R
# This preprocess shoud be used for mRNA, miRNA, lncRNA, CNV, methylation datasets
mRNA=read.csv("nsclc_mRNA_expr.csv",header=T, sep=",",row.names=1,check.names=F)
###Data filter
mRNA1=data.filter(mRNA,percentage=0.8)
####Data filling
mRNA1[mRNA1==0]<-NA
mRNA1<-as.matrix(mRNA1)
mRNA2=data.imputation(mRNA1,fun="mean")
###Feature Selection
mRNA3=SbyVar(mRNA2, cut.type="topk",value=1000)
###Data Standardization
mRNA4=data.normalization(mRNA3,type="feature_Median",log2=TRUE)
write.csv(mRNA4,file = "preprocessed_nsclc_mRNA.csv")
 ```
Code snippet: feature selection based on metabolism dataset
```
data(mRNAexp)
data_lipid=metabolism(mRNAexp,type="Amino_acid")
```
Code snippet:  Execute WMRCA based on five datasets
```R
data(mRNAexp)
data(miRNAexp)
data(lncRNAexp)
data(methylation)
data(CNVmatrix)
mRNAexp=data.filter(mRNAexp,percentage=0.6)
miRNAexp=data.filter(miRNAexp,percentage=0.6)
lncRNAexp=data.filter(lncRNAexp,percentage=0.6)
methylation=data.filter(methylation,percentage=0.6)
CNVmatrix=data.filter(CNVmatrix,percentage=0.6)
mRNAexp[mRNAexp==0]<-NA
miRNAexp[miRNAexp==0]<-NA
lncRNAexp[lncRNAexp==0]<-NA
methylation[methylation==0]<-NA
CNVmatrix[CNVmatrix==0]<-NA
mRNAexp<-as.matrix(mRNAexp)
miRNAexp<-as.matrix(miRNAexp)
lncRNAexp<-as.matrix(lncRNAexp)
methylation<-as.matrix(methylation)
CNVmatrix<-as.matrix(CNVmatrix)
mRNAexp=data.imputation(mRNAexp,fun="mean")
miRNAexp=data.imputation(miRNAexp,fun="mean")
lncRNAexp=data.imputation(lncRNAexp,fun="mean")
methylation=data.imputation(methylation,fun="mean")
CNVmatrix=data.imputation(CNVmatrix,fun="mean")
results = WMRCA(mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix,maxK=6,
reps=50,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",
innerLinkage="complete",seed=1262118388.71279,plot="pdf")
```
Code snippet: example for run Similarity Network Fusion (SNF)
```R
# runSNF
data(mRNAexp)
data(miRNAexp)
data(lncRNAexp)
data(methylation)
data(CNVmatrix)
result=runSNF(mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix, clusterNum=3, K=20, alpha=0.5, t=20)
result$group
# similarly, runiCluster, runCC, runCNMF
```

Code snippet: Visualization

```R
# survival analysis
data(results)
data(meta)
maxK = 6
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
 M = results[[i]]$consensusMatrix
 Fn = ecdf(M[lower.tri(M)])
 PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
Surplot(results,meta,optK=3)
# differential expression analysis
data("mRNAexp")
mRNAexp=as.matrix(mRNAexp)
Normal_Data=as.matrix(mRNAexp[,sample(1:100,20)])
result=Dif.limma(Tumor_Data=mRNAexp,Normal_Data=Normal_Data,group=NULL,topk=NULL,RNAseq=FALSE)
# proportion of ambiguous clustering (PAC) plot 
data(dot_df)
PACplot(dot_df)
# drugsensitivity (in development)
data(exprData)
data(studyResponse)
drugsensitivity(testMatrix=exprData,
studyResponse=studyResponse,
bortIndex=bortIndex,
drug="Bosutinib",
tissueType = "all",
batchCorrect = "eb",
selection=1,
dataset = "cgp2014",colours=c('#e94753','#47a4e9'))
# calculate indices
setwd("/Users/LGJ-15/Manuscripts/CFC_20230304")
inteMatrix=read.csv("integrated_results_12-Dec-2023 19.30.csv",header=T, sep=",",row.names=1, check.names=F)
a <- indices(inteMatrix,results_CoC,stas_all=TRUE,plot = TRUE,clusterNum=6,vote=TRUE)
# Check distribution
pdf(file="checkDistribution.pdf",height = 4,
    width = 4)
data.checkDistribution(mRNAexp)
# plot ROC
require(ROCR)
setwd("/Users/LGJ-15/Manuscripts/nsclc/results/cluster_subtype_info")
combined_results=read.csv("/Users/LGJ-15/Manuscripts/CFC_20230304/nsclc_mRNA_results_by_ConsensusClusterPlus_31-Dec-2023 03.22.csv",header=T, sep=",",check.names=F)
print(confusion_matrix  <- table(combined_results$mRNA_cc, combined_results$sample))
tn <- confusion_matrix [1, 1]
tp <- confusion_matrix [2, 2]
fp <- confusion_matrix [2, 1]
fn <- confusion_matrix [1, 2]
# # only for cnmf, The correct predictions are all on the anti-diagonal line.
# tn <- confusion_matrix [2, 1]
# tp <- confusion_matrix [1, 2]
# fp <- confusion_matrix [1, 1]
# fn <- confusion_matrix [2, 2]
sensitivity = tp/(tp + fn)# Calculate sensitivity, or TPR, true rate, recall rate (Recall)
specificity = tn/(tn + fp)
ppp = tp/(tp+fp)
npp = tn/(tn+fn)
accuracy = (tp+tn)/(tp+tn+fp+fn)#calculate accuracy
n = 3 # valitated float number 
result <- paste("Sensitivity = ", round(sensitivity, n) ,    #n is the number of bits reserved
                "\nSpecificity = ", round(specificity, n),       # /n  line breaks
                "\nPositive Predictive Value = ", round(ppp, n),
                "\nNegative Predictive Value = ", round(npp, n),
                "\nAccuracy = ", round(accuracy, n), "\n", sep="")
cat(result)
# Plot confusion matrix 
pdf(file="confusion matrix_mRNA_CC.pdf",height = 4,
    width = 4)
plot(confusion_matrix,col=c("#FF0000","#0000FF"),
      main="Confusion Matrix",
      xlab="Predicted",
      ylab="Actual"
      )
# Plot confusion matrix 
fourfoldplot(confusion_matrix,color=c("green","red"),main = "Confusion Matrix")
dev.off()
# Plot ROC curve
library(ROCR) 
library(RColorBrewer)
pred <- prediction(combined_results$mRNA_cc, combined_results$sample)
performance(pred,'auc')@y.values #AUC value
perf <- performance(pred,'tpr','fpr') # for CNMF, need to change the position of 'tpr' and 'fpr'
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
set1 = c(brewer.pal(9,"Set1"), 
         brewer.pal(8, "Dark2"))
pdf(file="ROC_mRNA_CC.pdf",height = 4,
    width = 4)
plot(perf,
     xlim =c(0,1), ylim=c(0,1),col=set1[2], 
     main = "",
     lwd = 3,
     cex.main=1.3,
     cex.lab=1.2, 
     cex.axis=1.2,
     font=1.2,
)
#plot(pred,col=set1[2],add = TRUE,lwd = 3)
legend("bottomright", bty = "n" ,col=set1[2], 
       lty=1, paste(c("AUC value ="),round(auc, 3)),lwd=3)
graphics::abline(a = 0,  b = 1,lwd = 2,lty=2)
dev.off()
# plot t-SNE
set.seed(123456)
intersection_tumor <- as.data.frame(lapply(mRNAexp, as.numeric))
row.names(intersection_tumor) <- row.names(mRNAexp)
colnames(intersection_tumor) <- colnames(mRNAexp)
intersection_tumor = na.omit(intersection_tumor)
intersection_tumor <- as.matrix(intersection_tumor)
library(Rtsne)
tSNE_res<- Rtsne(t(intersection_tumor), dims= 2, perplexity= 10, verbose= F, max_iter= 500, check_duplicates= F) 
tsne<- data.frame(tSNE1 = tSNE_res[["Y"]][,1], tSNE2= tSNE_res[["Y"]][,2], cluster = annCol$results) 
ggplot(tsne, aes(x= tSNE1, y= tSNE2, color= cluster)) + 
  geom_point(size= 4.5, alpha= 0.5) + 
  stat_ellipse(level= 0.85, show.legend = F) + 
  theme(legend.position= "top")+
  theme_bw()
```

## Features


### WMRCA implemented and integrated the state-of-the-art clustering algorithms and multi-omics data modeling analysis functionalities with R. 
#### 1) Key Points

  * We have refined the cluster-of-clusters method, specifically the Weighted Majority Rule-based Cluster-of-clusters Approach (WMRCA), for the integration of multi-omics data to accurately predict tumor subtypes.

  * The WMRCA model uses ten evaluation metrics to assess clustering performance and applies a weighted majority voting rule to automatically determine the optimal number of clusters. 

  * The WMRCA+ model, which features genes related to lipid metabolism, demonstrates superior performance compared to other commonly used clustering algorithms.

## Publication
Guojun Liu, et al. WMRCA: A weighted majority rule-based ﻿cluster-of-clusters approach for cancer subtype prediction using multi-omics data. *Briefings in Bioinformatics*, submitted.

## Contact

The WMRCA package is developed by Inner Mongolia Key Laboratory of Functional Genomics and Bioinformatics, Inner Mongolia University of Science and Technology, Baotou, China

  * Guojun Liu <gjliu77@gmail.com> 

