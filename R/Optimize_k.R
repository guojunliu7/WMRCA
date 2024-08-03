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
