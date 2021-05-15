library(dplyr)
library(reshape2)
library(qqman)

args=commandArgs(T)

args=c('DESEQ2normalized_merge_matrix.txt','genotype_overlap_peak.txt','0.00001','20individual_SNPs_only.pos.GWAS','cluster.txt')

signal <- read.table(args[1],header = T)

genotype <- read.table(args[2],header = F)
colnames(genotype) <- c('chrom1','start1','end1','SNP_id',colnames(signal),'chrom2','start2','end2')
genotype$PEAK_ID <- paste(genotype$chrom2,genotype$start2,genotype$end2,sep = '-')


#SNP information
SNP <- genotype[,4:(4+ncol(signal))]
SNP <- unique(SNP)
rownames(SNP) <- SNP$SNP_id
genotype <- genotype[,c(1:4,ncol(genotype))]


#caculate corelation
corelation <- c()
P <- c()
for (j in 1:nrow(genotype)){
  x=as.matrix(SNP[genotype$SNP_id[j],-1])
  y=as.matrix(signal[genotype$PEAK_ID[j],])
  rownames(x) <- NULL
  rownames(y) <- NULL
  
  cortable <- data.frame(x=as.vector(x),y=as.vector(y))  
  corR <- cor(cortable,method = "pearson")
  corP <- cor.test(x,y,method = "pearson")
  corP <- corP$p.value
  corelation <- append(corelation,corR[1,2])
  P <- append(P,corP)
}

genotype$CORRELATION <- corelation
genotype$P <- P
genotype$BH <- p.adjust(genotype$P)

pdf('qqplot.pdf')
qq(genotype$P)
dev.off()

sig <- genotype[genotype$P< as.numeric(args[3]),]


#####SNP annotation for individual specific
write.table(sig,file = 'SNP_associated_with_epi.txt',row.names=F,sep = "\t",quote = F)

#####GWAS annotation for individual specific

GWAS <- read.table(args[4],header = F)
colnames(GWAS) <- c("chrom",'pos','SNP_id','P_LM',"traits")
GWAS_ANNO <- inner_join(sig,GWAS,by='SNP_id')

write.table(GWAS_ANNO,file = 'GWAS_SNP_associated_with_epi.txt',row.names=F,sep = "\t",quote = F)

#####SNP annotation for ancestry specific
ancestry <- read.table(args[5],header = T)
sig_ancestry_ANNO <-inner_join(sig,ancestry,by="PEAK_ID")
write.table(sig_ancestry_ANNO,file = 'ancestry-specific-SNP_associated_with_epi.txt',row.names=F,sep = "\t",quote = F)



#####GWAS annotation for ancestry specific
GWAS_ancestry_ANNO <- inner_join(GWAS_ANNO,ancestry,by="PEAK_ID")
write.table(GWAS_ancestry_ANNO,file = 'ancestry-specific-GWAS_associated_with_epi.txt',row.names=F,sep = "\t",quote = F)

