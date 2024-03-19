load('pancancer/pancancer_drawdata.Rdata')
drawdata[1:5,1:10]
#Look at the number of tumors and normal tissues each
table(dplyr::filter(drawdata,Type=='Tumor')$Cancer)
table(dplyr::filter(drawdata,Type=='Normal')$Cancer)
library(ggpubr)
p <- ggboxplot(drawdata, x = "Cancer", y = "SPP1",
               fill = "Type",legend=F,palette =c("#87CEFA","#FFAEB9"))+
  theme_bw()
p
p1<-p + stat_compare_means(aes(group = Type), label = "p.signif")
dir.create("OUT")
ggsave("OUT/allcancer.pdf",plot = p1, width = 14, height = 7)


#############
genes<-c("GBM" ,"LGG","LIHC","CESC","LUAD","COAD","LAML","BRCA","ESCA","SARC","KIRP","STAD","PRAD",
         "SKCM","UCEC","HNSC","KIRC","LUSC","THYM","THCA","MESO","READ","PAAD","OV","TGCT","PCPG",
         "UVM" , "UCS" ,"BLCA","KICH","ACC","CHOL","DLBC")

###Combine the high and low grouping with the count data
datalist<-list()
for (i in 1:33) {
  file = paste0("E:/yjs/bioinformation/cancer/KM/High_Low/mixallcancer_result/",genes[i],".csv")
  m<-read.csv(file,row.names = 1)
  datalist[[i]] <-m[,c("ID","Cancer","clus0KM","clus1KM","clus2KM","clus3KM","clus4KM","clus5KM","clus6KM","clus7KM","clus8KM","clus9KM","clus10KM")] 
}
N<-do.call(rbind,datalist)
load("E:/yjs/bioinformation/cancer/KM/finaldata.Rdata")
data$ID<-rownames(data)
allcount<-merge(data,N,by = c("ID","Cancer"))
save(allcount,file = "allcount.Rdata")

################
load("allcount.Rdata")
genes<-unique(allcount$Cancer)
dir.create("cancers_count")
for (i in 1:33) {
  a<-allcount[which(allcount$Cancer==genes[i]),]
  files<-paste0("cancers_count/",genes[i],".csv")
  write.csv(a,file = files)
}


##########
group<-paste0("clus",0:10,"KM")
####Differential analysis of different subgroups and different types of cancer yielded differential genes

#####################
#clus0-SKCM LGG
# read data
gene<-c( "SKCM", "LGG")
diff<-list()
for (i in 1:length(gene)) {
  #Read the matrix
  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)
  #Correspond to groups and samples
  conditions<-readCount[,c("ID","clus0KM")]
  #Strip null values
  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))
  #Stripped of clinical information
  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]
  #Name the row and column names, and change the .in symbolname to -
  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus0KM)
  
  #Perform the cycle
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  # Run the Wilcoxon rank-sum test for each gene
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  #Deposit in the list
  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}
save(diff,file = "outfiles/diff_clus0.Rdata")

library("org.Hs.eg.db")
library(ReactomePA)
library(clusterProfiler)

result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}

save(result,file = "outfiles/result_clus0.Rdata")
#####################
#clus1-SKCM LGG
gene<-c("SKCM", "LGG")
diff<-list()
for (i in 1:length(gene)) {
  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)
  conditions<-readCount[,c("ID","clus1KM")]
  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))
  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]
  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus1KM)
  

  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  # Run the Wilcoxon rank-sum test for each gene
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")

  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

  
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}
save(diff,file = "outfiles/diff_clus1.Rdata")


result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file="outfiles/result_clus1.Rdata")

#####################
#clus2-SKCM LGG


files<-paste0("cancers_count/","SKCM",".csv")
readCount<-read.csv(file=files, header = T, row.names = 1)

conditions<-readCount[,c("ID","clus2KM")]

conditions <- conditions[complete.cases(conditions),]
b<-match(group,colnames(readCount))

readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
readCount1[1:4,1:4]

colname<-colnames(readCount1)
rowname<-readCount1[,1]
readCount1<-t(readCount1[,-1])
colnames(readCount1)<-rowname
colname<-gsub('[.]', '-', colname)
rownames(readCount1)<-colname[-1]
conditions<-as.factor(conditions$clus2KM)

# edger TMM normalize
library(edgeR)
y <- DGEList(counts=readCount1,group=conditions)
##Remove rows conssitently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)


pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")


conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

# Output results based on FDR threshold

outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
fdrThres=0.05
#Deposit in the list
diff<-outRst[outRst$FDR<fdrThres,]
save(diff,file = "diff_clus2.Rdata")


result<-list()
diffres<-diff
diffres$SYMBOL<-rownames(diffres)
SYMBOL<-diffres$SYMBOL
gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
diffres<-merge(gene_ID,diffres,by = "SYMBOL")
diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
diffres<-diffres[,c("ENTREZID","rank")]
geneList<-diffres[,2]
names(geneList) = as.character(diffres[,1])
geneList = sort(geneList, decreasing = TRUE)
result<- gsePathway(geneList,
                     organism= "human",
                     minGSSize= 10,
                     maxGSSize= 500,
                     pvalueCutoff= 0.05, 
                     pAdjustMethod= "BH",
                     verbose= FALSE,
                     eps= 0)
save(result,file = "result_clus2.Rdata")

#####################
#clus3-THYM SKCM LGG
gene<-c("THYM","SKCM", "LGG")
diff<-list()
for (i in 1:length(gene)) {

  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)

  conditions<-readCount[,c("ID","clus3KM")]

  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))

  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]

  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus3KM)
  

  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}
save(diff,file = "outfiles/diff_clus3.Rdata")


result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file = "outfiles/result_clus3.Rdata")

#####################
#clus4- THYM SKCM LGG 
gene<-c("THYM" ,"SKCM", "LGG")
diff<-list()
for (i in 1:length(gene)) {

  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)

  conditions<-readCount[,c("ID","clus4KM")]

  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))
  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]

  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus4KM)
  

  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  # Run the Wilcoxon rank-sum test for each gene
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05

  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}

save(diff,file = "outfiles/diff_clus4.Rdata")



result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file = "outfiles/result_clus4.Rdata")


#####################
#clus5-PAAD SKCM LGG THYM
gene<-c("PAAD", "LGG", "SKCM", "THYM")
diff<-list()
for (i in 1:length(gene)) {
  
  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)
  
  conditions<-readCount[,c("ID","clus5KM")]
  
  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))
  
  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]
  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus5KM)
  
  
  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}
save(diff,file = "outfiles/diff_clus5.Rdata")


result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file = "outfiles/result_clus5.Rdata")


#####################
#clus6-SKCM LGG GBM 
gene<-c("SKCM", "LGG", "GBM")
diff<-list()
for (i in 1:length(gene)) {

  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)

  conditions<-readCount[,c("ID","clus6KM")]

  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))

  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]
 
  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus6KM)
  

  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
   
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05

  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}
save(diff,file = "diff_clus6.Rdata")



result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file = "outfiles/result_clus6.Rdata")


#####################
#clus7-UCEC BRCA SKCM LGG 
gene<-c("UCEC", "BRCA", "SKCM", "LGG")
diff<-list()
for (i in 1:length(gene)) {

  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)

  conditions<-readCount[,c("ID","clus7KM")]

  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))

  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]
  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus7KM)
  

  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
   
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}
save(diff,file = "outfiles/diff_clus7.Rdata")



result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file = "outfiles/result_clus7.Rdata")


#####################
#clus8-SKCM LGG
gene<-c( "SKCM" ,"LGG")
diff<-list()
for (i in 1:length(gene)) {

  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)

  conditions<-readCount[,c("ID","clus8KM")]

  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))

  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]

  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus8KM)
  

  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}
save(diff,file = "outfiles/diff_clus8.Rdata")



result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file = "outfiles/result_clus8.Rdata")


#####################
#clus9-LUAD PAAD KIRP THYM MESO LGG ACC UCEC
gene<-c("LUAD" ,"PAAD", "KIRP", "THYM", "MESO", "LGG", "ACC","UCEC")
diff<-list()
for (i in 1:length(gene)) {
 
  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)
  
  conditions<-readCount[,c("ID","clus9KM")]
 
  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))
  
  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]
  
  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus9KM)
  
  
  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  
  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}

save(diff,file = "outfiles/diff_clus9.Rdata")


result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file = "outfiles/result_clus9.Rdata")


#####################
#clus10-SKCM LGG ACC UVM
gene<-c("SKCM", "LGG")
diff<-list()
for (i in 1:length(gene)) {
  
  files<-paste0("cancers_count/",gene[i],".csv")
  readCount<-read.csv(file=files, header = T, row.names = 1)
  
  conditions<-readCount[,c("ID","clus10KM")]
  
  conditions <- conditions[complete.cases(conditions),]
  b<-match(group,colnames(readCount))
  
  readCount1<-readCount[match(conditions$ID,readCount$ID),-c(b,2:10)]
  readCount1[1:4,1:4]
  colname<-colnames(readCount1)
  rowname<-readCount1[,1]
  readCount1<-t(readCount1[,-1])
  colnames(readCount1)<-rowname
  colname<-gsub('[.]', '-', colname)
  rownames(readCount1)<-colname[-1]
  conditions<-as.factor(conditions$clus1
  # edger TMM normalize
  library(edgeR)
  y <- DGEList(counts=readCount1,group=conditions)
  ##Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  ##Perform TMM normalization and transfer to CPM (Counts Per Million)
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  
  outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  fdrThres=0.05
  
  diff[[i]]<-outRst[outRst$FDR<fdrThres,]
}
save(diff,file = "outfiles/diff_clus10.Rdata")



result<-list()
for (i in 1:length(diff)) {
  diffres<-diff[[i]]
  diffres$SYMBOL<-rownames(diffres)
  SYMBOL<-diffres$SYMBOL
  gene_ID <- bitr(SYMBOL, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  diffres<-merge(gene_ID,diffres,by = "SYMBOL")
  diffres$rank = (-log(diffres$pValues))*diffres$log2foldChange
  diffres<-diffres[,c("ENTREZID","rank")]
  geneList<-diffres[,2]
  names(geneList) = as.character(diffres[,1])
  geneList = sort(geneList, decreasing = TRUE)
  ret_ges<- gsePathway(geneList,
                       organism= "human",
                       minGSSize= 10,
                       maxGSSize= 500,
                       pvalueCutoff= 0.05, 
                       pAdjustMethod= "BH",
                       verbose= FALSE,
                       eps= 0)
  result[[i]]<-ret_ges
}
save(result,file = "outfiles/result_clus10.Rdata")




#clus0
rm(list = ls())
load("outfiles/result_clus0.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result
cancer5<-result[[5]]@result
cancer6<-result[[6]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4,cancer5,cancer6)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))


library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
#Select the up-down pathway of interest, and select the up-down pathway below
df = allcancer
df = df[order(df$enrichmentScore),]
#Set up groups：1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p0<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p0,filename = "OUT/clus0.pdf",height = 15,width =20)


############
##clus1
rm(list = ls())
load("outfiles/result_clus1.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4)

library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))

df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p1<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p1,filename = "OUT/clus1.pdf",height = 15,width =20)

############
#clus2
rm(list = ls())
load("outfiles/result_clus2.Rdata")
cancer1<-result@result

library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)

df = cancer1
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p2<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p2,filename = "OUT/clus2.pdf",height = 15,width =20)



############
#clus3
rm(list = ls())
load("outfiles/result_clus3.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result

cancers<-rbind(cancer1,cancer2,cancer3)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))


library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p3<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p3,filename = "OUT/clus3.pdf",height = 15,width =20)


############
#clus4
rm(list = ls())
load("outfiles/result_clus4.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))


library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p4<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p4,filename = "OUT/clus4.pdf",height = 15,width =20)


############
#clus5
rm(list = ls())
load("outfiles/result_clus5.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))


library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p5<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p5,filename = "OUT/clus5.pdf",height = 15,width =20)

############
#clus6
rm(list = ls())
load("outfiles/result_clus6.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))


library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p6<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p6,filename = "OUT/clus6.pdf",height = 15,width =20)


############
#clus7
rm(list = ls())
load("outfiles/result_clus7.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result
cancer5<-result[[5]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4,cancer5)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))


library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p7<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p7,filename = "OUT/clus7.pdf",height = 10,width =20)

############
#clus8
rm(list = ls())
load("outfiles/result_clus8.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))


library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p8<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p8,filename = "OUT/clus8.pdf",height = 15,width =20)

############
#clus9
rm(list = ls())
load("outfiles/result_clus9.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result
cancer5<-result[[5]]@result
cancer6<-result[[6]]@result
cancer7<-result[[7]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4,cancer5,cancer6,cancer7)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))


library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p9<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p9,filename = "OUT/clus9.pdf",height = 15,width =20)

############
#clus10
rm(list = ls())
load("outfiles/result_clus10.Rdata")
cancer1<-result[[1]]@result
cancer2<-result[[2]]@result
cancer3<-result[[3]]@result
cancer4<-result[[4]]@result

cancers<-rbind(cancer1,cancer2,cancer3,cancer4)

library(tidyverse)
allcancer<-cancers %>% group_by(Description)  %>% slice(which.min(qvalues))

library(ggplot2)
library(forcats)
library(ggstance)
library(ggstatsplot)
df = allcancer
df = df[order(df$enrichmentScore),]
#1--ES>0;   -1--ES<0
up = head(subset(df, enrichmentScore>0),10);up$group=1
down = tail(subset(df, enrichmentScore<0),10);down$group=-1
dat=rbind(up,down)
#dat$group = factor(dat$group);str(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
p10<-ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity", aes(fill=factor(group, levels = c(1,-1), labels = c("ES>0","ES<0")))) + 
  xlab("Pathway names") +
  ylab("-log10(P.adj)") +
  coord_flip() + 
  theme_ggstatsplot() +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4),
                     labels=c("4", "2", "0","2","4")) +
  scale_fill_manual(values = c("tomato","#0068D2")) + 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle("Pathway Enrichment") 
ggsave(p10,filename = "OUT/clus10.pdf",height = 15,width =20)




##Pathway enrichment analysis of PD-L1 on characteristic genes 
library(readxl)
matrix<-read.csv("GSE181815/GSE181815_20180514_ThymusCa_countMatrix_anno.csv",check.names = F)
#Grouping information
pheno <- read_excel("GSE181815/GSE181815_ClinicalData.xlsx")
pheno<-pheno[,c(1,8)]
colnames(pheno)<-c("sample","group")

pheno[which(pheno$group == "Complete Response"),2] <- "Response"
pheno[which(pheno$group == "Partial Response"),2] <- "Response"
pheno[which(pheno$group == "Progressive Disease"),2] <- "No-Response"
pheno[which(pheno$group == "Stable Disease"),2] <- "No-Response"
##Name the expression matrix row names
matrix<-matrix[,-1]
matrix <- aggregate(x = matrix[,-10],
                    by = list(matrix$GeneSymbol),
                    FUN = max)
rownames(matrix)<-matrix$Group.1
matrix<-matrix[,-1]
colnames(matrix)<-gsub('.genes.results', '', colnames(matrix))
colnames(matrix)<-gsub('[.]', '-', colnames(matrix))
pheno<-pheno[-10,]#Sample excluded from analysis due to poor quality
matrix<-matrix[,pheno$sample]
#Remove cells with full expression of 0
matrix <- matrix[rowMeans(matrix)>0,]   
group_list<-pheno$group



library(edgeR)

dgelist <- DGEList(counts = matrix, group = group_list)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')


design <- model.matrix(~group_list)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)

fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(lrt, 'GSE181815/control_treat.glmLRT.txt', sep = '\t', col.names = NA, quote = FALSE)


fit <- glmQLFit(dge, design, robust = TRUE)
lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
write.table(lrt, 'GSE181815/control_treat.glmQLFit.txt', sep = '\t', col.names = NA, quote = FALSE)


gene_diff <- read.delim('GSE181815/control_treat.glmLRT.txt', row.names = 1, sep = '\t', check.names = FALSE)
#First, sort the table in ascending order by FDR value, and continue to sort by log2FC descending under the same FDR value
gene_diff <- gene_diff[order(gene_diff$FDR, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]

#log2FC≥1 & FDR<0.01  up
#log2FC≤-1 & FDR<0.01 down
#The rest identify none, which represents non-differential genes
gene_diff[which(gene_diff$logFC >= 1 & gene_diff$FDR < 0.01),'sig'] <- 'up'
gene_diff[which(gene_diff$logFC <= -1 & gene_diff$FDR < 0.01),'sig'] <- 'down'
gene_diff[which(abs(gene_diff$logFC) <= 1 | gene_diff$FDR >= 0.01),'sig'] <- 'none'

#Output a summary table of differential genes for selection
diff<-gene_diff[which(gene_diff$PValue<0.05),]
save(diff,file = "GSE181815/diff.Rdata")


library(dplyr)
library(clusterProfiler)
load("GSEA_PATH.Rdata")
diff <- diff %>% arrange(desc(logFC))
geneList<-diff$logFC#Extract the foldchange from largest to smallest
names(geneList) <- rownames(diff)#Add the corresponding Symbol to the foldchange extracted above

gsea <- GSEA(geneList,
             minGSSize = 1,
             pvalueCutoff = 1,
             seed = 123,
             TERM2GENE = path) 

save(gsea,file = "GSE181815/gsea.Rdata")
head(gsea)
a<-data.frame(gsea)
a<-a[a$pvalue<0.05,]
a<-a[c(2,6,7,3,4,5,1),]
library(pheatmap)
library(enrichplot)
library(ggplot2)
genes<-a$core_enrichment
genes<-strsplit(a$core_enrichment, "/")

group<-paste0("clus",c(2:4,6:8,10))
for (i in 1:7) {
  p1<-gseaplot2(gsea,
                group[i], 
                color="red", 
                title = group[i],
                base_size = 10, 
                subplots = 1:3, 
                pvalue_table = F) 
  filename1 = paste0("OUT/GSE181815/gsea",group[i],".png")
  ggsave(p1,filename = filename1,width = 8,height = 8)
  name<-match(genes[[i]],rownames(diff))
  b<-diff[name,]
  c<-data.frame(b$logFC)
  rownames(c)<-rownames(b)
  p2<-pheatmap(c,show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F,border_color = NA,color = colorRampPalette(c("white", "red"))(1000))
  filename2 = paste0("OUT/GSE181815/heatmap",group[i] ,".png")
  ggsave(p2,filename =filename2 ,width = 2,height = 8)
}




#################
library(data.table)
matrix<-fread("GSE126044/GSE126044_counts.txt.gz",check.names = F)

pheno <- read.table("GSE126044/SraRunTable.txt",header = T,sep = ",")
pheno<-pheno[,c('Sample.Name','patient_response')]
pheno1<-data.frame(names=c("Dis_01","Dis_02","Dis_11","Dis_12","Dis_15","Dis_04","Dis_16","Dis_17","Dis_03","Dis_10","Dis_06","Dis_18","Dis_09","Dis_05","Dis_07","Dis_08"),
          Sample.Name=c("GSM3589673","GSM3589674","GSM3589675","GSM3589676","GSM3589677","GSM3589678","GSM3589679","GSM3589680","GSM3589681","GSM3589682","GSM3589683","GSM3589684","GSM3589685","GSM3589686","GSM3589687","GSM3589688"))
pheno <- merge(pheno,pheno1,by = "Sample.Name")
pheno<-pheno[,c(2,3)]
colnames(pheno)<-c("group","sample")
pheno<-pheno[order(pheno$group),]

dim(matrix)
#[1] 18747    17
matrix <- aggregate(x = matrix[,-1],
                    by = list(matrix$V1),
                    FUN = max)#
rownames(matrix)<-matrix$Group.1
matrix<-matrix[,-1]
matrix<-matrix[,pheno$sample]
matrix <- matrix[rowMeans(matrix)>0,]   
group_list<-pheno$group


###
library(edgeR)

dgelist <- DGEList(counts = matrix, group = group_list)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

design <- model.matrix(~group_list)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(lrt, 'GSE126044/control_treat.glmLRT.txt', sep = '\t', col.names = NA, quote = FALSE)

fit <- glmQLFit(dge, design, robust = TRUE)
lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
write.table(lrt, 'GSE126044/control_treat.glmQLFit.txt', sep = '\t', col.names = NA, quote = FALSE)

gene_diff <- read.delim('GSE126044/control_treat.glmLRT.txt', row.names = 1, sep = '\t', check.names = FALSE)
gene_diff <- gene_diff[order(gene_diff$FDR, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]

#log2FC≥1 & FDR<0.01  up
#log2FC≤-1 & FDR<0.01  down
#The rest identify none, which represents non-differential genes
gene_diff[which(gene_diff$logFC >= 1 & gene_diff$FDR < 0.01),'sig'] <- 'up'
gene_diff[which(gene_diff$logFC <= -1 & gene_diff$FDR < 0.01),'sig'] <- 'down'
gene_diff[which(abs(gene_diff$logFC) <= 1 | gene_diff$FDR >= 0.01),'sig'] <- 'none'

diff<-gene_diff[gene_diff$PValue<0.05,]
save(diff,file = "GSE126044/diff.Rdata")

library(dplyr)
library(clusterProfiler)
load("GSEA_PATH.Rdata")
diff <- diff %>% 
  arrange(desc(logFC))
geneList<-diff$logFC
names(geneList) <- rownames(diff)

gsea <- GSEA(geneList,
             minGSSize = 1,
             seed = 123,
             pvalueCutoff = 1,
             TERM2GENE = path) 

save(gsea,file = "GSE126044/gsea.Rata")
head(gsea)
a<-data.frame(gsea)
a<-a[a$p.adjust<0.05,]
a<-a[c(1,2,3,9,10,6,4,8,5,7),]
library(pheatmap)
library(enrichplot)
library(ggplot2)
genes<-a$core_enrichment
genes<-strsplit(a$core_enrichment, "/")


group<-paste0("clus",c(0:8,10))
for (i in 1:10) {
  p1<-gseaplot2(gsea,
                group[i], 
                color="red",
                title = group[i],
                base_size = 10, 
                subplots = 1:3, 
                pvalue_table = F) 
  filename1 = paste0("OUT/GSE126044/gsea",group[i],".png")
  ggsave(p1,filename = filename1,width = 8,height = 8)
  name<-match(genes[[i]],rownames(diff))
  b<-diff[name,]
  c<-data.frame(b$logFC)
  rownames(c)<-rownames(b)
  p2<-pheatmap(c,show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F,border_color = NA,color = colorRampPalette(c("white", "red"))(1000))
  filename2 = paste0("OUT/GSE126044/heatmap",group[i] ,".png")
  ggsave(p2,filename =filename2 ,width = 2,height = 16)
}




#################

library(readxl)
matrix<-read_excel("GSE78220/GSE78220_PatientFPKM.xlsx")
pheno <- read.table("GSE78220/SraRunTable.txt",header = T,sep = ",")
pheno<-pheno[,c('anti.pd.1_response','GEO_Accession..exp.')]
colnames(pheno)<-c("group","sample")
pheno[which(pheno$group == "Complete Response"),1] <- "Response"
pheno[which(pheno$group == "Partial Response"),1] <- "Response"
pheno[which(pheno$group == "Progressive Disease"),1] <- "No-Response"
号
sample_ID <- read_excel("GSE78220/sample_ID.xlsx")
colnames(matrix)<-gsub(".baseline","",colnames(matrix))
colnames(matrix)<-gsub(".OnTx","",colnames(matrix))
match(colnames(matrix),sample_ID$id)
pheno<-merge(pheno,sample_ID,by = "sample")
pheno<-pheno[,c(2,3)]

dim(matrix)
#[1] 25268    29
matrix <- aggregate(x = matrix[,-1],
                    by = list(matrix$Gene),
                    FUN = max)#
rownames(matrix)<-matrix$Group.1
matrix<-matrix[,-1]
pheno<-pheno[order(pheno$group),]
n<-match(pheno$id,colnames(matrix))
matrix<-matrix[,n]
matrix <- matrix[rowMeans(matrix)>0,]   
group_list<-pheno$group


###
expMatrix <- matrix
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)

tpms[1:3,]
colSums(tpms)
## Enforce the qualifying order
group_list <- factor(group_list,levels = c("No-Response","Response"),ordered = F)
exprSet <- tpms
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#Determine if the data needs to be transformed
exprSet <- log2(exprSet+1)
#Variance Analysis:
dat <- exprSet
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 
deg$FDR<-p.adjust(deg$P.Value,method="fdr",n=length(deg$P.Value))
#save(deg,file = 'deg.Rdata')
diff<-deg[which(deg$P.Value<0.05),]
save(diff,file = "GSE78220/diff.Rdata")

library(dplyr)
library(clusterProfiler)
load("GSEA_PATH.Rdata")
diff <- diff %>% arrange(desc(logFC))
geneList<-diff$logFC
names(geneList) <- rownames(diff)

gsea <- GSEA(geneList,
             minGSSize = 1,
             seed = 123,
             pvalueCutoff = 1,
             TERM2GENE = path) 

save(gsea,file = "GSE78220/gsea.Rdata")
head(gsea)
a<-data.frame(gsea)
a<-a[which(a$p.adjust<0.05),]
a<-a[c(2,1,6,4,3,5,7),]
library(pheatmap)
library(enrichplot)
library(ggplot2)
genes<-a$core_enrichment
genes<-strsplit(a$core_enrichment, "/")

group<-paste0("clus",c(0,1,3,5,6,8,10))
for (i in c(1:7)) {
  p1<-gseaplot2(gsea,
                group[i], 
                color="red", 
                title =  group[i],
                base_size = 10, 
                subplots = 1:3,
                pvalue_table = F) 
  filename1 = paste0("OUT/GSE78220/gsea", group[i],".png")
  ggsave(p1,filename = filename1,width = 8,height = 8)
  name<-match(genes[[i]],rownames(diff))
  b<-diff[name,]
  c<-data.frame(b$logFC)
  rownames(c)<-rownames(b)
  p2<-pheatmap(c,show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F,border_color = NA,color = colorRampPalette(c("white", "red"))(1000))
  filename2 = paste0("OUT/GSE78220/heatmap", group[i] ,".png")
  ggsave(p2,filename =filename2 ,width = 2,height = 8)
}

