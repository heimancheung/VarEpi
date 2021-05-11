library(DESeq2)
library(dplyr)


##########################################################
################individual-specific peak##################
##########################################################


args=commandArgs(T)
#args <- c('D:/ExperimentData/20+2rice_epi/pipeline/H3K4me3_seedlings+ZH11_split_ancestry.txt','10','3',"0.001","15")

readcounttable <- read.table(args[1],header = T)


#transform dataframe
row.names(readcounttable) <- readcounttable[,1]
readcounttable <- readcounttable[,-1]


#get factor from readcounttable col name
splitList <- strsplit(colnames(readcounttable),split = '_')
repVetor <- c()
for (i in c(1:ncol(readcounttable))){
  repVetor[i] <- paste(splitList[[i]][1],splitList[[i]][2],sep="_") 
}
rm(splitList)

#creat deseq2 matrix
condition <- factor(repVetor,levels = repVetor[2*(1:(ncol(readcounttable)))])

#get deseq2 normalized matrix
colData <- data.frame(row.names=colnames(readcounttable), condition)
dds <- DESeqDataSetFromMatrix(readcounttable, colData, design= ~ condition)
dds <- estimateSizeFactors(dds)
dds2 <- counts(dds,normalized=T)




#calculate anova 
##open an empty dataframe with rownames and colnames
ANOresData <- data.frame(matrix(NA,nrow(dds2),3))  
rownames(ANOresData) <- rownames(dds2)
colnames(ANOresData) <-  c("bP","F","P")

for (i in 1:nrow(dds2)){
  ANO_matrix <- data.frame(X=as.vector(dds2[i,]),A=factor(repVetor, levels = repVetor[2*1:18]))
  bres <- bartlett.test(X~A,ANO_matrix)
  aovres <- summary.aov(aov(X~A,data = ANO_matrix))
  ANOresData[i,1] <- bres$p.value
  ANOresData[i,2] <- aovres[[1]]$`F value`[1]
  ANOresData[i,3] <- aovres[[1]]$`Pr(>F)`[1]
}

variable <- ANOresData[ANOresData$F> as.integer(args[2]),]
conserved <- ANOresData[ANOresData$F<=as.integer(args[2]),]

write.table(variable,file="variable_ano_res.txt",quote = F,sep='\t')
write.table(conserved,file="non_variable_ano_res.txt",quote = F,sep='\t')


####save data set
splitList <- strsplit(colnames(dds2),split = '_')
repVetor <- c()
for (i in 1:(ncol(dds2))){
  repVetor[(i)] <- paste(splitList[[i]][1], splitList[[i]][2],sep = "_")
}
rm(splitList)

mergedds2 <- dds2[,seq(2,ncol(dds2),2)] + dds2[,seq(1,ncol(dds2),2)]
colnames(mergedds2) <- repVetor[seq(1,length(repVetor),2)]

write.table(dds2,file='DESEQ2normalized_matrix.txt',quote = F,sep='\t')
write.table(mergedds2,file='DESEQ2normalized_merge_matrix.txt',quote = F,sep='\t')



########################################################
################ancestry-specific peak##################
########################################################


####ancestry-specific peak
library("preprocessCore")
require(lattice)
library(dplyr)
library(isva)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(factoextra)
library(preprocessCore)



######merge data
data <- readcounttable 
  
splitList <- strsplit(colnames(data),split = '_')
repVetor <- c()
for (i in 1:(ncol(data))){
  repVetor[(i)] <- paste(splitList[[i]][1], splitList[[i]][2],sep = "_")
}
rm(splitList)

mergeData <- data[,seq(2,ncol(data),2)] + data[,seq(1,ncol(data),2)]
colnames(mergeData) <- repVetor[seq(1,length(repVetor),2)]
rownames(mergeData) <- rownames(data)

####normaized raw signal
asinhdata <- asinh(mergeData)
QNdata <- normalize.quantiles(as.matrix(asinhdata),copy=TRUE)
QNdata <- as.data.frame(QNdata)
colnames(QNdata) <- colnames(mergeData)
rownames(QNdata) <- rownames(mergeData)

######variable peak CV top 60%
QNdata_var <- QNdata[rownames(variable),]

CV <- apply(t(QNdata_var), 2,sd) / apply(t(QNdata_var),2,mean)
data_CV <- QNdata_var[order(CV,decreasing = T),]
data_CV <- data_CV[1:floor(nrow(data_CV)*0.6),]



####To indentify group-specific regions, ISVA removing confounding factor
splitList <- strsplit(colnames(data_CV),split = '_')
ancestry <- c()
for (i in 1:ncol(data_CV)){
  ancestry[i] <- splitList[[i]][1] 
}
rm(splitList)

pheno.v <- ancestry

data.m <- data_CV

prinR <- princomp(data.m)
pdf('./screeplot.pdf')
screeplot(prinR)
dev.off()

ls = list()
for (RANDOM in 1:20){
  ncomp = as.integer(args[3])
  
  lm.o <- lm(t(data.m) ~ as.factor(pheno.v))
  res.m <- t(lm.o$res)											
  model <- model.matrix(~1 + as.factor(pheno.v))
  
  fICA.o <- fastICA(res.m, n.comp = ncomp)
  tmp.m <- t(fICA.o$A)									#A: estimated mixing matrix
  isv.m <- tmp.m
  sd <- 1/sqrt(ncol(data.m) - 3)							#??????????????????
  
  for (k in 1:ncol(tmp.m)) {
    cor.v <- as.vector(cor(t(data.m), tmp.m[, k]))		
    z.v <- 0.5 * log((1 + cor.v)/(1 - cor.v))
    pv.v <- 2 * pnorm(abs(z.v), 0, sd, lower.tail = FALSE)
    tmp.s <- sort(pv.v, decreasing = FALSE, index.return = TRUE)
    qv.o <- qvalue(pv.v)								
    nsig <- length(which(qv.o$qvalues < 0.05))
    if (nsig < 500) {
      nsig <- 500
    }
    red.m <- data.m[tmp.s$ix[1:nsig], ]				
    fICA.o <- fastICA(red.m, n.comp = ncomp)
    cor.v <- abs(cor(tmp.m[, k], t(fICA.o$A)))	
    kmax <- which.max(cor.v)
    isv.m[, k] <- t(fICA.o$A)[, kmax]
    print(paste("Built ISV ", k, sep = ""))
  }
  
  
  
  isva.o <- isv.m
  
  selisv.idx <- 1:ncol(isva.o)
  pv.m <- NULL
  
  selisv.m <- matrix(isva.o[, selisv.idx], ncol = length(selisv.idx))
  mod <- model.matrix(~pheno.v + selisv.m)
  modNULL <- model.matrix(~selisv.m)
  df1 <- dim(mod)[2]
  df0 <- dim(modNULL)[2]
  pv.v <- rep(0, nrow(data.m))
  Id <- diag(ncol(data.m))
  #resid <- data.m %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  lm1 <- lm(t(data.m) ~ as.factor(pheno.v) + selisv.m)
  resid2 <- lm1$residuals
  
  rss1 <- rowSums(t(resid2) * t(resid2))
  rm(resid2)
  #residNULL <- data.m %*% (Id - modNULL %*% solve(t(modNULL) %*% 
  #                                                 modNULL) %*% t(modNULL))
  lm2 <- lm(t(data.m) ~  selisv.m)
  residNULL <- lm2$residuals
  
  rssNULL <- rowSums(t(residNULL) * t(residNULL))
  rm(residNULL)
  fstats <- ((rssNULL - rss1)/(df1 - df0))/(rss1/(ncol(data.m) - 
                                                    df1))
  pv.v <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (ncol(data.m) - 
                                                     df1))
  pv.s <- sort(pv.v, decreasing = FALSE, index.return = TRUE)
  qv.v <- qvalue(pv.s$x)$qvalue
  ntop <- length(which(qv.v < as.numeric(args[4])))
  
  
  pred.idx <- pv.s$ix[1:ntop]
  lm.o <- lm(t(data.m[pred.idx, ]) ~ pheno.v + selisv.m)
  tstats.v <- unlist(lapply(summary(lm.o), function(x) {
    x$coeff[2, 3]
  }))
  
  
  lm.m <- cbind(tstats.v, pv.s$x[1:ntop], qv.v[1:ntop])
  colnames(lm.m) <- c("t-stat", "P-value", "q-value")
  ls[[RANDOM]] <- lm.m
}

interLS <- rownames(ls[[1]])
for (i in 2:10){
  interLS <- intersect(interLS,rownames(ls[[i]]))
}

interLS <- gsub('[Response ]', '', interLS)

write.table(interLS,file = 'ancestry_specific.peak',
           quote = F,sep = '\t',row.names = F)



###############################################
################visualization##################
###############################################

SVA <- data.m[interLS,]
SVA$PEAK_ID <- rownames(SVA)
SVA[,-ncol(SVA)] <- normalize.quantiles(as.matrix(SVA[,-ncol(SVA)]),copy=TRUE)

rowM <- apply(SVA[,-ncol(SVA)],1,mean)
rowSD <- apply(SVA[,-ncol(SVA)],1,sd)
SVA[,-ncol(SVA)] <- (SVA[,-ncol(SVA)]-rowM) / rowSD
SVA <- na.omit(SVA)


###group number 
WSS <- c()
for (k in 1:50){
  km_result <- kmeans(SVA[,-ncol(SVA)],k)
  SSE = km_result$tot.withinss
  times = 1
  while(times <= 10){
    km_result <- kmeans(SVA[,-ncol(SVA)],k)
    if (km_result$tot.withinss < SSE) {SSE <- km_result$tot.withinss}
    times <- times + 1
  }
  WSS <- append(WSS,SSE)
}

pdf('Within_groups_sum_of_squares.pdf')
plot(1:50,WSS,xlab="number of cluster",ylab = "Within groups sum of squares")
dev.off()


###kmeans clstering
km_result <- kmeans(SVA[,-ncol(SVA)],as.numeric(args[5]))


km_clu <- as.data.frame(km_result$cluster)
km_clu$PEAK_ID <- rownames(km_clu)

KMCluM <- inner_join(SVA,km_clu,by='PEAK_ID')
KMCluM <- KMCluM[order(KMCluM$`km_result$cluster`),]
row.names(KMCluM) <- KMCluM$PEAK_ID
KMCluM <- KMCluM[,-ncol(KMCluM)]

##ordering
clu.m <- melt(KMCluM)

xOr <- colnames(KMCluM)
xOr <- xOr[-length(xOr)]
xOr <- xOr[order(xOr)]

clu.m$PEAK_ID <- factor(rep(KMCluM$PEAK_ID,each=1,time=length(xOr)),levels = KMCluM$PEAK_ID)
clu.m$variable <-  factor(clu.m$variable,levels = xOr)


###plotting
pdf('./ancestry_specific_peak_heatmap.pdf')
ggplot(clu.m) + geom_tile(aes(x=variable,y=PEAK_ID,fill=value))  +
  scale_fill_gradient2(low = '#1565c0',mid ='#fffde7' ,high='#d50000',midpoint = 0.5) +
  theme(panel.background=element_blank())  + xlab('rice') + ylab('region') +
  theme_grey(base_size=8) + theme(axis.text.y = element_blank()) 
dev.off()

colnames(km_clu) <- c('cluster','PEAK_ID')

write.table(km_clu[,c(2,1)],'cluster.txt',quote = F,sep = '\t',row.names = F)
