indices <- function(inteMatrix,results_CoC,stas_all = TRUE , plot = TRUE,clusterNum=6,vote=TRUE) {
  if (stas_all == TRUE){
    requireNamespace("ade4", quietly = TRUE)
    requireNamespace("fpc", quietly = TRUE)
    requireNamespace("cluster", quietly = TRUE)
    if (!dir.exists("output")){
      dir.create("output")
    }else{
      print("dir exists")
    }
    maxK= clusterNum # important
    inteMatrix_square_matrix <- dist(t(inteMatrix),method ="euclidean") # # Compute pairwise-distance matrices, generate a square matrix, cluster.stats functions requirs square matrix
    # hc_stats <- cluster.stats(mRNAexp1_stats,results_CoC[[i]]$consensusClass)
    # hc_stats
    stcl = lapply(2:maxK, function(i) cluster.stats(inteMatrix_square_matrix,results_CoC[[i]]$consensusClass) )  # set 
    ldunn = sapply(1:(maxK-1), function(i) stcl[[i]]$dunn )
    ldunn2 = sapply(1:(maxK-1),function(i) stcl[[i]]$dunn2)
    lwbr  = sapply(1:(maxK-1), function(i) stcl[[i]]$wb.ratio ) #c(stcl.B$wb.ratio,stcl.hc$wb.ratio,stcl.Whc$wb.ratio,stcl.W2hc$wb.ratio,stcl.km$wb.ratio)
    lch = sapply(1:(maxK-1), function(i) stcl[[i]]$ch ) #c(stcl.B$ch,stcl.hc$ch,stcl.Whc$ch,stcl.W2hc$ch,stcl.km$ch)
    lentr = sapply(1:(maxK-1), function(i) stcl[[i]]$entropy )
    lpear = sapply(1:(maxK-1),function(i) stcl[[i]]$pearsongamma)
    lgap = sapply(1:(maxK-1),function(i) stcl[[i]]$widestgap)
    lsil = sapply(1:(maxK-1),function(i) stcl[[i]]$avg.silwidth)
    ## Mean cophenetic distance
    cdl = lapply(2:maxK, function(i) as.dist(1-results_CoC[[i]]$consensusMatrix) )
    corl =sapply(cdl, cor, inteMatrix_square_matrix)
    ## Script name: PCA 
    maxK = length(results_CoC)  #
    Kvec = 2:maxK
    x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
    PAC = rep(NA,length(Kvec))
    names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
    for(i in Kvec){
      M = results_CoC[[i]]$consensusMatrix
      Fn = ecdf(M[lower.tri(M)])   #提取出共识矩阵下三角的数据，然后将用ecdf 方法生成拟合曲线
      PAC[i-1] = Fn(x2) - Fn(x1)   #计算0.1到0.9之间的面积
    }
    optK = Kvec[which.min(PAC)]    #面积最小值对应K为最佳K
    cat("The optimal K is",optK,"according to PAC\n")
    cat("The optimal K is",which.max(lsil)+1,"according to Mean Silhouette\n")
    cat("The optimal K is",which.max(lch)+1,"according to Calinski-Harabasz index\n")
    cat("The optimal K is",which.max(lentr)+1,"according to Entropy\n")
    cat("The optimal K is",which.max(ldunn)+1,"according to Dunn\n")
    cat("The optimal K is",which.max(ldunn2)+1,"according to Dunn 2\n")
    cat("The optimal K is",which.min(lwbr)+1,"according to Within-between SS ratio\n")
    cat("The optimal K is",which.max(corl)+1,"according to Mean cophenetic distance\n")
    cat("The optimal K is",which.max(lpear)+1,"according to Pearson gamma\n")
    cat("The optimal K is",which.max(lgap)+1,"according to Gap\n")
    #################################################################
    # Script name: plot Silhouette with function "silhouette"
    #################################################################
    # pdf(paste("output/Silhouette",".pdf",sep=""),h=8,w=4*(maxK-1))
    # par(mfrow=c(1,(maxK-1)))
    # lsil2 = vector("list",(maxK-1))
    # for(i in 2:maxK){
    #   sil = silhouette(results[[i]]$consensusClass,dist(t(mRNAexp1),method = "euclidean"))
    #   sizes = table(results[[i]]$consensusClass)
    #   plot(sil,col=rep(results[[i]]$clrs[[3]],rep=sizes) ,main=paste("K=",i))
    #   lsil2[[i-1]]=sil
    # }
    # dev.off()
    # msil = sapply(1:(maxK-1), function(i) mean( lsil2[[i]][,3] ) )
    # cdl = lapply(2:maxK, function(i) as.dist(1-results[[i]]$consensusMatrix ) )
    # md = dist(t(mRNAexp1),method = "euclidean")
    # corl =sapply(cdl, cor, md)
    ###################################################################
    ## Script name: # plot clustering stats
    ###################################################################
    pdf(paste("output/Cluster_separation_stats",".pdf",sep=""),h=3.5,w=3.5*(maxK-1))
    par(mfrow=c(1,(maxK-1)),family="Times")
    co = rep(1,(maxK-1))
    co[which.min(PAC)]=2
    barplot(PAC, names.arg = paste("K =",2:maxK), las=2, main = "Proportion of ambiguous clustering",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3)
    co = rep(1,(maxK-1))
    co[which.max(lsil)]=2
    barplot(lsil,names.arg = paste("K =",2:maxK) ,las=2,main = "Mean Silhouette",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3)
    co = rep(1,(maxK-1))
    co[which.max(lch)]=2
    barplot(lch ,names.arg = paste("K =",2:maxK) ,las=2,main = "Calinski-Harabasz index",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3)
    co = rep(1,(maxK-1))
    co[which.max(lentr)]=2
    barplot(lentr, names.arg = paste("K =",2:maxK), las=2, main = "Entropy",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3)
    co = rep(1,(maxK-1))
    co[which.max(ldunn)]=2
    barplot(ldunn ,names.arg = paste("K =",2:maxK) ,las=2,main = "Dunn index",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3)
    co = rep(1,(maxK-1))
    co[which.max(ldunn2)]=2
    barplot(ldunn2 ,names.arg = paste("K =",2:maxK) ,las=2,main = "Dunn index 2",col= co,density = 50,cex.names = 1.30,cex.axis = 1.3)
    co = rep(1,(maxK-1))
    co[which.min(lwbr)]=2
    barplot(lwbr ,names.arg = paste("K =",2:maxK) ,las=2,main = "Within-between SS ratio",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3)
    co = rep(1,(maxK-1))
    co[which.max(corl)]=2
    barplot(corl,names.arg = paste("K =",2:maxK) ,las=2,main = "Mean cophenetic distance",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3)
    co = rep(1,(maxK-1))
    co[which.max(lpear)]=2
    barplot(lpear,names.arg = paste("K =",2:maxK) ,las=2,main = "Pearson gamma",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3,cex.lab=1.5)
    co = rep(1,(maxK-1))
    co[which.max(lgap)]=2
    barplot(lgap,names.arg = paste("K =",2:maxK) ,las=2,main = "Gap",col= co,density = 50,cex.names = 1.3,cex.axis = 1.3,cex.lab=1.5)
    # mtext("Gap", side=2, line=2.2, cex=2)
    dev.off()
    
    if (vote==TRUE){
      ###################################################################
      ## Script name: Majority voting
      ###################################################################
      pdf(paste("output/Majority_voting",".pdf",sep=""),h=4,w=3.6)
      # par(mfrow=c(1,(maxK-2)),family="Times")
      Optimal_k_list <- c(optK,
                          which.max(lsil)+1,
                          which.max(lch)+1,
                          which.max(lentr)+1,
                          which.max(ldunn)+1,
                          which.max(ldunn2)+1,
                          which.min(lwbr)+1,
                          which.max(corl)+1,
                          which.max(lpear)+1,
                          which.max(lgap)+1)

      names(which.max(table(Optimal_k_list)))
      requireNamespace("RColorBrewer", quietly = TRUE)
      coul <- brewer.pal(5, "Set2") 
      barplot(table(Optimal_k_list),col=coul,ylab="Majority voting",cex.axis = 1.3,cex.names = 1.3,cex.lab=1.5)
      dev.off()
    }
    Optimal_k <- names(which.max(table(Optimal_k_list)))
    cat("....................................\n")
    cat("The final optimal K list is ",Optimal_k_list," according to Majority rule\n")
    cat("The final optimal K is",Optimal_k,"according to Majority rule\n")
  }
}
optimize <- function(method,plot=FALSE){
  if (method==CPCC) {
    d2=""
    cor1=""
    cor <- numeric()
    for(i in 2:length(results)){
      k2 <- results[[i]]$consensusMatrix
      d1 <- dist(k2)
      hc <- hclust(d1, "complete")
      d2 <- cophenetic(hc)
      cor <- c(cor, assign(paste0("cor", i-1), cor(d1, d2)))
    }
    k <- numeric()
    for(i in 2:length(results)) {
      k <- c(k,k=i) 
    }
    if (plot==TRUE) {
      pdf("CPCC.pdf")
      plot(k, cor,
           # col ="red", 
           # pch=22,
           # xlim=c(0,6),
           # ylim=c(0,5),
           # lwd=2,
           # xlab="WEE",
           ylab = "Cophenetic Correlation Coefficient",
           type="b"
      )
      dev.off()     
    }
  }
  return(cor)
}

# 
# pdf("PAC_CoC.pdf")
# maxK = length(results)  #
# Kvec = 2:maxK
# x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
# PAC = rep(NA,length(Kvec)) 
# names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
# for(i in Kvec){
#   M = results[[i]]$consensusMatrix
#   Fn = ecdf(M[lower.tri(M)])   #提取出共识矩阵下三角的数据，然后将用ecdf 方法生成拟合曲线
#   PAC[i-1] = Fn(x2) - Fn(x1)   #计算0.1到0.9之间的面积
# }
# plot(PAC)
# dev.off()
# optK = Kvec[which.min(PAC)]    #面积最小值对应K为最佳K
# optK
# PAC