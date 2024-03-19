###################
##First, select the characteristic genes for each cluster
load("cancer/macrophage/Files/scRNA_acc_markers.Rdata")
DefaultAssay(object = scRNA) <- "integrated"

c0 <- FindMarkers(scRNA, ident.1=0, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c1 <- FindMarkers(scRNA, ident.1=1, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c2 <- FindMarkers(scRNA, ident.1=2, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c3 <- FindMarkers(scRNA, ident.1=3, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c4 <- FindMarkers(scRNA, ident.1=4, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c5 <- FindMarkers(scRNA, ident.1=5, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c6 <- FindMarkers(scRNA, ident.1=6, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c7 <- FindMarkers(scRNA, ident.1=7, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c8 <- FindMarkers(scRNA, ident.1=8, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c9 <- FindMarkers(scRNA, ident.1=9, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
c10 <- FindMarkers(scRNA, ident.1=10, slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")

library(dplyr)
c0 = c0 %>% subset(c0$p_val_adj<0.05 & abs(c0$avg_log2FC) > 0.3)
c1 = c1 %>% subset(c1$p_val_adj<0.05 & abs(c1$avg_log2FC) > 0.3)
c2 = c2 %>% subset(c2$p_val_adj<0.05 & abs(c2$avg_log2FC) > 0.3)
c3 = c3 %>% subset(c3$p_val_adj<0.05 & abs(c3$avg_log2FC) > 0.3)
c4 = c4 %>% subset(c4$p_val_adj<0.05 & abs(c4$avg_log2FC) > 0.3)
c5 = c5 %>% subset(c5$p_val_adj<0.05 & abs(c5$avg_log2FC) > 0.3)
c6 = c6 %>% subset(c6$p_val_adj<0.05 & abs(c6$avg_log2FC) > 0.3)
c7 = c7 %>% subset(c7$p_val_adj<0.05 & abs(c7$avg_log2FC) > 0.3)
c8 = c8 %>% subset(c8$p_val_adj<0.05 & abs(c8$avg_log2FC) > 0.3)
c9 = c9 %>% subset(c9$p_val_adj<0.05 & abs(c9$avg_log2FC) > 0.3)
c10 = c10 %>% subset(c10$p_val_adj<0.05 & abs(c10$avg_log2FC) > 0.3)
#dir.create("Findmarker")

write.csv(c0,file = "Findmarker/cluster0marker.csv")
write.csv(c1,file = "Findmarker/cluster1marker.csv")
write.csv(c2,file = "Findmarker/cluster2marker.csv")
write.csv(c3,file = "Findmarker/cluster3marker.csv")
write.csv(c4,file = "Findmarker/cluster4marker.csv")
write.csv(c5,file = "Findmarker/cluster5marker.csv")
write.csv(c6,file = "Findmarker/cluster6marker.csv")
write.csv(c7,file = "Findmarker/cluster7marker.csv")
write.csv(c8,file = "Findmarker/cluster8marker.csv")
write.csv(c9,file = "Findmarker/cluster9marker.csv")
write.csv(c10,file = "Findmarker/cluster10marker.csv")




##################
## KM    ANOVA  ##

library(TCGAbiolinks)
projects <- getGDCprojects()
library(dplyr)
projects <- projects %>% 
  as.data.frame() %>% 
  dplyr::select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))

gdcdata=function(i){
  ## 0.Operational information
  print(paste0("Downloading number ",i,",project name: ",projects$project_id[i]))
  ## 1.Query information
  query.exp = GDCquery(project = projects$project_id[i], 
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")
  ## 2. download
  GDCdownload(query.exp)
  ## 3.Merge multiple data
  pre.exp = GDCprepare(query = query.exp)
  ## 4.Extract expression data
  library(SummarizedExperiment)
  countsdata = SummarizedExperiment::assay(pre.exp,1)
  fpkmdata=SummarizedExperiment::assay(pre.exp,5)
  tpmdata=SummarizedExperiment::assay(pre.exp,4)
  gene_id=data.frame(id=rowData(pre.exp)@listData[["gene_id"]], gene_name= rowData(pre.exp)@listData[["gene_name"]],gene_type=rowData(pre.exp)@listData[["gene_type"]])
  counts=cbind(gene_id,countsdata)
  fpkm=cbind(gene_id,fpkmdata)
  tpm=cbind(gene_id,tpmdata)
  #Clinical information
  clinical <- GDCquery_clinic(project = projects$project_id[i], type = "clinical")
  ## 5.Save data
  filename1 = paste0("result/",projects$project_id[i],"-counts.txt")
  filename2 = paste0("result/",projects$project_id[i],"-fpkm.txt")
  filename3 = paste0("result/",projects$project_id[i],"-tpm.txt")
  filename4 = paste0("result/",projects$project_id[i],"-clinical.txt")
  write.table(counts,filename1,sep="\t",col.names=T,row.names=F,quote=F) 
  write.table(fpkm,filename2,sep="\t",col.names=T,row.names=F,quote=F) 
  write.table(tpm,filename3,sep="\t",col.names=T,row.names=F,quote=F) 
  write.table(clinical,filename4,sep="\t",col.names=T,row.names=F,quote=F) 
}
dir.create("result")
for (i in 1:3) {
  gdcdata(i)
}







####Read clinical data
projects <- getGDCprojects()
projects <- projects %>% 
  as.data.frame() %>% 
  select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))
clin<-list()
clin1<-list()
for (i in 1:33) {
  filename4 = paste0("result/",projects$project_id[i],"-clinical.txt")
  a <- read.table(filename4,header = T,sep = "\t")
  a <-a %>% 
    dplyr::select(bcr_patient_barcode, 
                  vital_status,
                  days_to_death, days_to_last_follow_up,
                  gender,age_at_index)
  a$Cancer = projects$tumor[i]
  
  #1、Calculate time to live
  a$days_to_death[is.na(a$days_to_death)] <- 0   #NA are marked as 0
  a$days_to_last_follow_up[is.na(a$days_to_last_follow_up)] <- 0
  a$days=apply(cbind(a$days_to_death, a$days_to_last_follow_up), 1, max)
  #Time is recorded in months, with two decimal places
  a$time=round(a$days/30,2)
  
  #2、According to the definition of life and death living is 0 and dead is 1
  a$event=ifelse(a$vital_status=='Alive',0,1)
  clin1[[i]]<-a
}

clin2<-do.call("rbind",clin1)
save(clin2,file = "clindata.Rdata")


#####
source("http://bioconductor.org/biocLite.R")
biocLite("maftools")

library(maftools)

exprs<-list()
library(tidyr)
for (i in 1:33) {
  filename1 = paste0("result/",projects$project_id[i],"-counts.txt")
  mRNA<-read.table(filename1,header = T,check.names = F)
  mRNA$Cancer = projects$tumor[i]
  mRNA<-mRNA[,-c(1,3)]
  mRNA[,2:ncol(mRNA)]<-lapply(mRNA[,2:ncol(mRNA)],as.factor)#Converting to a factor and then a value does not force NA to be added
  mRNA[,2:ncol(mRNA)]<-lapply(mRNA[,2:ncol(mRNA)],as.numeric)
  data <- aggregate(x = mRNA[,-1],
                    by = list(mRNA$gene_name),
                    FUN = mean)
  name<-substr(colnames(data), 1, 12)
  data<-rbind(name,data)
  exprs[[i]]<-data

}
save(exprs,file = "exprs.Rdata")
exprs1<-list()
for (i in 1:33) {
  a<-exprs[[i]]
  a<-a[,-which(colnames(a)=="Cancer")]
  tumor<-rep(projects$tumor[i],ncol(a))
  a<-rbind(tumor,a)
  exprs1[[i]]<-a
}

mRNA<-do.call("cbind",exprs1)
save(mRNA,file = "mRNA.Rdata")

mRNA[1:4,1:7]
mRNA[1,1]<-"Cancer"
a<-which(colnames(mRNA)=='Group.1')
a<-a[-1]
mRNA<-mRNA[,-a]
rownames(mRNA)<-mRNA[,1]
mRNA<-mRNA[,-1]
mRNA<-mRNA[-c(3:5),]#Remove non-coding rRNA
mRNA1<-rbind(colnames(mRNA),mRNA)
save(mRNA1,file = "mRNA1.Rdata")

mRNA<-data.frame(t(mRNA1))

group_list=ifelse(as.numeric(substr(rownames(mRNA),14,15)) < 10,'tumor','normal')
group_list<-data.frame(group_list)
group_list$ID<-rownames(mRNA)
head(group_list)

#Consolidation of data
mRNA$ID<-row.names(mRNA)
mRNA=inner_join(mRNA,group_list,by ="ID",copy=T)
mRNA[1:5,1:5]
mRNA<-select(mRNA,ID,Group.1,Cancer,group_list,everything())
mRNA[1:5,1:5]
save(mRNA,file = "counts.Rdata")

#####There are uncleaned row names for each tumor data
a<-which(mRNA$A1CF=="A1CF")
mRNA[a[1:4],1:9]
a<-a[-1]
count<-mRNA[-a,]
colnames(count)[5:59428]<-mRNA[175,5:59428]
save(count,file = "count.Rdata")



load("clindata.Rdata")
count[1:3,1:6]
colnames(count)[2]<-"bcr_patient_barcode"
count<-count[which(count$group_list=="tumor"),]
data<-merge(clin2,count,by="bcr_patient_barcode")
data<-data[,-c(1,3,4)]
data[1:3,1:9]
rownames(data)<-data$ID
data<-data[,-9]
colnames(data)[1]<-"status"
colnames(data)[3]<-"Age"
colnames(data)[4]<-"Cancer"

save(data,file = "finaldata.Rdata")




#######
load("pancancer/gtf.Rdata")
load("finaldata.Rdata")
count<-data[,-c(1:8)]
count[1:4,1:5]
count<-cbind(rownames(count),count)
count1<-data.frame(t(count))
colnames(count1)<-count1[1,]
count1[1:4,1:5]
count1<-count1[-1,]
save(count1,file="linshi.Rdata")


#Select Start and End to calculate the gene length to obtain a .tpm file
a=dplyr::filter(gtf,type=="gene",gene_biotype=="protein_coding") 
b=dplyr::select(a,c(gene_name,start,end))#The three selected in the gtf are mainly the required gene ID and the data replaced by the ensemble ID
length<-abs(b$end-b$start)
c<-cbind(length,b)
c<-c[,1:2]
gene_name<-rownames(count1)
mRNA<-cbind(gene_name,count1)
d=merge(c,mRNA,by ="gene_name")#One-to-one mapping by genename
if(!require("dplyr")) BiocManager::install("dplyr")
library("dplyr")

#TPM
mRNAdata=distinct(d,gene_name,.keep_all = T)
eff_length <- mRNAdata[,"length"]/1000
count<-mRNAdata[,c(-1,-2)]
rownames(count)<-mRNAdata[,1]
count[1:5,1:5]
count[,2:ncol(mRNA)]<-lapply(mRNA[,2:ncol(mRNA)],as.factor)
count<-data.frame(lapply(count,as.numeric))
rownames(count)<-mRNAdata[,1]
a<-colnames(mRNAdata)
a<-a[-(1:2)]
colnames(count)<-a
x <- count / eff_length
mat_tpm <- t( t(x) / colSums(x) ) * 1e6 
rownames(mat_tpm)<- rownames(count) 
mat_tpm[1:4,1:4]
save(mat_tpm,file = "TPM.Rdata")



load("TPM.Rdata")
#score
#In order to ensure that the random numbers are consistent, here is the seed sequence
set.seed(123)
#Calculate the mean for each row
rowmean=apply(mat_tpm,1,mean)
#Calculate the standard deviation for each row
rowsd=apply(mat_tpm,1,sd)
data1=sweep(mat_tpm,1,rowmean)
data2=sweep(data1,1,rowsd,'/')
mRNA<-data.frame(t(data2))
mRNA[1:3,1:6]
load("finaldata.Rdata")
data[1:3,1:10]
ID<-rownames(data)
data<-data[,c(1:8)]
data<-cbind(ID,data)
ID<-rownames(mRNA)
mRNA<-cbind(ID,mRNA)
zdata<-merge(data,mRNA,by = "ID")


c0<-read.csv("Findmarker/cluster0marker.csv",row.names = 1)
c1<-read.csv("Findmarker/cluster1marker.csv",row.names = 1)
c2<-read.csv("Findmarker/cluster2marker.csv",row.names = 1)
c3<-read.csv("Findmarker/cluster3marker.csv",row.names = 1)
c4<-read.csv("Findmarker/cluster4marker.csv",row.names = 1)
c5<-read.csv("Findmarker/cluster5marker.csv",row.names = 1)
c6<-read.csv("Findmarker/cluster6marker.csv",row.names = 1)
c7<-read.csv("Findmarker/cluster7marker.csv",row.names = 1)
c8<-read.csv("Findmarker/cluster8marker.csv",row.names = 1)
c9<-read.csv("Findmarker/cluster9marker.csv",row.names = 1)
c10<-read.csv("Findmarker/cluster10marker.csv",row.names = 1)

######
#Characteristic genes are averaged
clus0<-rownames(c0)
a<-vector()
for (i in 1:length(clus0)) {
  a[i]<-match(clus0[i],colnames(zdata))
}
zdata$clus0<-rowMeans(select(zdata,na.omit(a)))
clus0<-colnames(zdata)[na.omit(a)]
clus1<-rownames(c1)
a<-vector()
for (i in 1:length(clus1)) {
  a[i]<-match(clus1[i],colnames(zdata))
}
zdata$clus1<-rowMeans(select(zdata,na.omit(a)))
clus1<-colnames(zdata)[na.omit(a)]
clus2<-rownames(c2)
a<-vector()
for (i in 1:length(clus2)) {
  a[i]<-match(clus2[i],colnames(zdata))
}
zdata$clus2<-rowMeans(select(zdata,na.omit(a)))
clus2<-colnames(zdata)[na.omit(a)]
clus3<-rownames(c3)
a<-vector()
for (i in 1:length(clus3)) {
  a[i]<-match(clus3[i],colnames(zdata))
}
clus3<-colnames(zdata)[na.omit(a)]
zdata$clus3<-rowMeans(select(zdata,na.omit(a)))
clus4<-rownames(c4)
a<-vector()
for (i in 1:length(clus4)) {
  a[i]<-match(clus4[i],colnames(zdata))
}
zdata$clus4<-rowMeans(select(zdata,na.omit(a)))
clus4<-colnames(zdata)[na.omit(a)]
clus5<-rownames(c5)
a<-vector()
for (i in 1:length(clus5)) {
  a[i]<-match(clus5[i],colnames(zdata))
}
zdata$clus5<-rowMeans(select(zdata,na.omit(a)))
clus5<-colnames(zdata)[na.omit(a)]
clus6<-rownames(c6)
a<-vector()
for (i in 1:length(clus6)) {
  a[i]<-match(clus6[i],colnames(zdata))
}
zdata$clus6<-rowMeans(select(zdata,na.omit(a)))
clus6<-colnames(zdata)[na.omit(a)]
clus7<-rownames(c7)
a<-vector()
for (i in 1:length(clus7)) {
  a[i]<-match(clus7[i],colnames(zdata))
}
zdata$clus7<-rowMeans(select(zdata,na.omit(a)))
clus7<-colnames(zdata)[na.omit(a)]
clus8<-rownames(c8)
a<-vector()
for (i in 1:length(clus8)) {
  a[i]<-match(clus8[i],colnames(zdata))
}
zdata$clus8<-rowMeans(select(zdata,na.omit(a)))
clus8<-colnames(zdata)[na.omit(a)]
clus9<-rownames(c9)
a<-vector()
for (i in 1:length(clus9)) {
  a[i]<-match(clus9[i],colnames(zdata))
}
zdata$clus9<-rowMeans(select(zdata,na.omit(a)))
clus9<-colnames(zdata)[na.omit(a)]
clus10<-rownames(c10)
a<-vector()
for (i in 1:length(clus10)) {
  a[i]<-match(clus10[i],colnames(zdata))
}
zdata$clus10<-rowMeans(select(zdata,na.omit(a)))
clus10<-colnames(zdata)[na.omit(a)]
dataTemp <- list(clus0,clus1,clus2,clus3,clus4,clus5,clus6,clus7,clus8,clus9,clus10)
dataTemp2 <- data.frame(do.call(cbind, 
                                lapply(lapply(dataTemp, unlist), `length<-`, 
                                       max(lengths(dataTemp)))))
colnames(dataTemp2)<-paste0("clus",0:10)
write.csv(dataTemp2,file = "Findmarker/signaturegenes_2.csv")


save(zdata,file = "zdata_clus.Rdata")





load("zdata_clus.Rdata")
load("reclindata.Rdata")
clind$ID<-rownames(clind)
data<-merge(clind,zdata,by=c("ID","Cancer"))
zdata<-data[,-c(8,9)]
colnames(zdata)[3]<-"days"
colnames(zdata)[4]<-"time"

#See the expression of individual genes in pan-cancer
#The gene with the highest differential expression was obtained
library(ggplot2)
p0<-ggplot(zdata,aes(x = reorder(Cancer, clus0, FUN = median), y = clus0,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p1<-ggplot(zdata,aes(x = reorder(Cancer, clus1, FUN = median), y = clus1,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p2<-ggplot(zdata,aes(x = reorder(Cancer, clus2, FUN = median), y = clus2,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p3<-ggplot(zdata,aes(x = reorder(Cancer, clus3, FUN = median), y = clus3,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p4<-ggplot(zdata,aes(x = reorder(Cancer, clus4, FUN = median), y = clus4,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p5<-ggplot(zdata,aes(x = reorder(Cancer, clus5, FUN = median), y = clus5,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p6<-ggplot(zdata,aes(x = reorder(Cancer, clus6, FUN = median), y = clus6,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p7<-ggplot(zdata,aes(x = reorder(Cancer, clus7, FUN = median), y = clus7,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p8<-ggplot(zdata,aes(x = reorder(Cancer, clus8, FUN = median), y = clus8,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p9<-ggplot(zdata,aes(x = reorder(Cancer, clus9, FUN = median), y = clus9,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
p10<-ggplot(zdata,aes(x = reorder(Cancer, clus10, FUN = median), y = clus10,color=Cancer))+
  geom_boxplot()+
  theme_bw()+labs(x="Cancer")+
  theme(legend.position='none',axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
#dir.create("OUT")
#dir.create("OUT/KM")
#dir.create("OUT/Exprs")

ggsave("OUT/Exprs/mixallcancer_result/cluster0.pdf",plot = p0, width = 14, height = 7)


rm(list = ls())

library(survival)
library(survminer)
load("zdata_clus.Rdata")
zdata[which(is.na(zdata$days)),c(3,4)]<-0
zdata[1:4,1:10]
genes<-unique(zdata$Cancer)
genes
gene<-c('clus0','clus1','clus2','clus3','clus4','clus5','clus6','clus7','clus8','clus9','clus10')
#gene <- 'SPP1'
##The genes of each tumor were grouped and calculated, and the matrix was plotted with dot plots
#########
#HR
genes
for (p in 1:length(genes)) {
  m<-zdata[which(zdata$Cancer==genes[p]),]
  a<-matrix(nrow = nrow(m),ncol = length(gene))
  #Cycle through the grouping
  for (i in 1:length(gene)) {
    for (j in 1:nrow(m)) {
      if (m[j,gene[i]]<= quantile(as.numeric(m[,which(colnames(m)==gene[i])]))[[2]]) {
        a[j,i] = 'Low'
      }
      else if(m[j,gene[i]] > quantile(as.numeric(m[,which(colnames(m)==gene[i])]))[[4]]){
        a[j,i] = 'High'
      }
    }
  }
  
  a<-as.data.frame(a)
  colnames(a)<-paste0(gene,"KM")
  #Combine data and grouping information
  b<-cbind(m,a)
  files <- paste0("High_Low/mixallcancer_result/",genes[p])
  write.csv(b,file = paste0(files,".csv"))
  #res.cox <- coxph(Surv(time, status) ~ JCHAIN, data = b)
  #summary(res.cox)
  ##All into one table
  covariates <- gene
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(days, event)~', x)))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = b)})
  # Extract the data and make a data table
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res<-as.data.frame(res)
  name <- paste0("cox_cancer/mixallcancer_result/",genes[p])
  write.csv(res,file = paste0(name,".csv"))
}



#############
#Merge all files
load("zdata_clus.Rdata")
zdata[1:4,1:10]
genes<-unique(zdata$Cancer)
genes
gene<-c('clus0','clus1','clus2','clus3','clus4','clus5','clus6','clus7','clus8','clus9','clus10')

genes
listdata<-list()
for (i in 1:length(genes)) {
  files<- paste0("cox_cancer/mixallcancer_result/",genes[i])
  a<-read.csv( paste0(files,".csv"),header = T)
  a$cancer <- genes[i]
  a<-a[,c(1,3,5,6)]
  colnames(a)<-c("cluster","HR","P_value","Cancer")
  a$HR<-gsub("\\(.*\\)","",a$HR)
  listdata[[i]]<-a
}

for (i in 1:(length(listdata)-1)) {
  listdata[[i+1]]<-rbind(listdata[[i]],listdata[[i+1]])
}
data<-as.data.frame(listdata[[33]])
which(is.na(data$P_value))
head(data)
data$FDR<-p.adjust(data$P_value,method="fdr",length(data$P_value))


#It is divided into two matrices: HR and FDR
library(tidyr)
data1<-data[,c(1,2,4)]
data1<-spread(data1,cluster, HR)
rownames(data1)<-data1$Cancer


data2<-data[,c(1,4,5)]
data2<-spread(data2,cluster, FDR)
rownames(data2)<-data2$Cancer


data1[,-1]<-as.data.frame(lapply(data1[,-1],as.numeric))
data2[,-1]<-as.data.frame(lapply(data2[,-1],as.numeric))
data1<-data1[,-1]
data2<-data2[,-1]
head(data1)
head(data2)

#All FDR values are judged, and those with p<0.01 are marked with "**", and those with p-value 0.01<p<0.05 are marked with "*".
if (!is.null(data2)){
  ssmt <- data2< 0.01
  data2[ssmt] <-'**'
  smt <- data2 >0.01& data2 <0.05
  data2[smt] <- '*'
  data2[!ssmt&!smt]<- ''
} else {
  data2 <- F
}


data1<-data.frame(t(data1))
data1<-data1[c(1:2,4:11,3),]
data2<-data.frame(t(data2))
data2<-data2[c(1:2,4:11,3),]

bk =  c(seq(0,0.9,by=0.01),seq(1,10,by = 0.1))
mycol<-colorRampPalette(c("#0068D2","white","tomato"))(length(bk)-1)
#When HR is 1, data1 is 0
library(pheatmap)
pheatmap(data1, scale = "none",border=NA,breaks = bk,
         cluster_rows = F,cluster_cols = T,
         display_numbers = data2,fontsize_number = 12,
         number_color = "black",color=mycol,
         angle_col = 0,fontsize_col=10,
         height = 7,width = 17,filename = "E:\\yjs\\bioinformation\\cancer\\Plots\\p3\\cox1.png")

pheatmap(data1, scale = "none",border=NA,breaks = bk,
         cluster_rows = F,cluster_cols = T,
         legend_breaks = c(0, 1, 10),legend_labels = c("0", "1", "10"),
         display_numbers = data2,fontsize_number = 12,
         number_color = "black",color=mycol,
         angle_col = 0,fontsize_col=10,
         height = 7,width = 17,filename = "E:\\yjs\\bioinformation\\cancer\\Plots\\p3\\cox2.png")

#############
#Survival analysis
load("zdata_clus.Rdata")
zdata[1:4,1:10]
genes<-unique(zdata$Cancer)
genes
gene<-c('clus0','clus1','clus2','clus3','clus4','clus5','clus6','clus7','clus8','clus9','clus10')


library(survival)
library(survminer)
genekm<-paste0(gene,"KM")
listdata<-list()
p<-matrix(nrow = 1,ncol = 33)
for (i in 1:length(genes)) {
  files<- paste0("High_Low/SPP1/",genes[i])
  a<-read.csv( paste0(files,".csv"),header = T,row.names = 1)
  listdata[[i]]<-a
  ### Batch drawing
  for (j in 1:length(genekm)){
    group <- a[[genekm[j]]]
    diff = survdiff(Surv(days , event )~ group ,data = a)
    pValue = 1-pchisq(diff$chisq,df=1)
    p[j,i]<-pValue
    fit <- survfit(Surv(days, event) ~ group,data = a)
    surPlot = ggsurvplot(fit,
                         data = a,
                         conf.int = T,
                         legend.labs = c("High","Low"),
                         legend.title = genekm[j],
                         xlab = "Time",
                         ylab = "Overall survival",
                         #break.time.by = 1,
                         risk.table.title = "",
                         palette = c("tomato","#0068D2"),
                         ggtheme =theme_light(),
                         risk.table = T,
                         risk.table.height = .25)
    names<-paste(genes[i],genekm[j],sep = "_")
    png(file = paste("E:\\yjs\\bioinformation\\cancer\\Plots\\SPP1\\",names,".png",sep = ),
        width = 620, 
        height = 620,          
        units = "px",          
        bg = "white",          
        res = 300)              
    print(surPlot)
    dev.off()
    
  }
}
####FDR is calculated for the characteristic gene p value of each tumor
p<-data.frame(p)
colnames(p)<-genes
FDR<-matrix(nrow=11,ncol=33)
for (i in 1:ncol(p)) {
  FDR[,i]<-p.adjust(p[,i],method="fdr",length(p[,i]))
}
colnames(FDR)<-genes
write.csv(p,file = "OUT/surv/mixallcancer_result/pvalue.csv")
write.csv(FDR,file = "OUT/surv/mixallcancer_result/FDR.csv")




##############
load('zdata_clus.Rdata')
zdata[1:4,1:10]
genes<-unique(zdata$Cancer)
genes
gene<-c('clus0','clus1','clus2','clus3','clus4','clus5','clus6','clus7','clus8','clus9','clus10')
library(tableone)

####GBM
b<-read.csv("High_Low/mixallcancer_result/GBM.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus6
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus6 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/GBM_clus6.csv")

####LUAD
b<-read.csv("High_Low/mixallcancer_result/LUAD.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LUAD_clus9.csv")

####BRCA
b<-read.csv("High_Low/mixallcancer_result/BRCA.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus7
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus7 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/BRCA_clus7.csv")

####ACC
b<-read.csv("High_Low/mixallcancer_result/ACC.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/ACC_clus9.csv")

####KIRP
b<-read.csv("High_Low/mixallcancer_result/KIRP.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/KIRP_clus9.csv")

####LIHC
b<-read.csv("High_Low/mixallcancer_result/LIHC.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LIHC_clus9.csv")

####MESO
b<-read.csv("High_Low/mixallcancer_result/MESO.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/MESO_clus9.csv")

####UCEC
b<-read.csv("High_Low/mixallcancer_result/UCEC.csv",header = T,row.names = 1)
b[1:3,1:8]
#clus7
res.cox <- coxph(Surv(days, event) ~ Age + clus7 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/UCEC_clus7.csv")

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/UCEC_clus9.csv")

####PAAD
b<-read.csv("High_Low/mixallcancer_result/PAAD.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus5
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus5 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/PAAD_clus5.csv")

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/PAAD_clus9.csv")

####THYM
b<-read.csv("High_Low/mixallcancer_result/THYM.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus3
res.cox <- coxph(Surv(days, event) ~ Age + gender+ clus3 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/THYM_clus3.csv")

#clus4
res.cox <- coxph(Surv(days, event) ~ Age + gender+ clus4 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/THYM_clus4.csv")

#clus5
res.cox <- coxph(Surv(days, event) ~ Age + gender+ clus5 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/THYM_clus5.csv")

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + gender+ clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/THYM_clus9.csv")

####SKCM
b<-read.csv("High_Low/mixallcancer_result/SKCM.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus0
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus0, data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus0.csv")

#clus1
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus1 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus1.csv")

#clus2
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus2 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus2.csv")

#clus3
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus3 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus3.csv")

#clus4
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus4 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus4.csv")

#clus5
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus5 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus56.csv")

#clus6
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus6 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus6.csv")

#clus7
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus7 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus7.csv")

#clus8
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus8 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus8.csv")

#clus10
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus10 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/SKCM_clus10.csv")


###LGG
b<-read.csv("High_Low/mixallcancer_result/LGG.csv",header = T,row.names = 1)
b[1:3,1:8]

#clus0
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus0 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus0.csv")

#clus1
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus1 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus1.csv")

#clus2
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus2 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus2.csv")

#clus3
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus3 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus3.csv")

#clus4
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus4 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus4.csv")

#clus5
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus5 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus5.csv")

#clus6
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus6 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus6.csv")

#clus7
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus7 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus7.csv")

#clus8
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus8 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus8.csv")

#clus9
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus9 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus9.csv")

#clus10
res.cox <- coxph(Surv(days, event) ~ Age + gender+clus10 , data = b)
res.cox1<-summary(res.cox)
colnames(res.cox1$conf.int)
multi1<-as.data.frame(round(res.cox1$conf.int[, c(1, 3, 4)], 2))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(res.cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Merge the results of the two extractions into a table, and name it result
result <-cbind(multi1,multi2);result
#The row name is changed to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
write.csv(result,file = "COXdata/mixallcancer_result/LGG_clus10.csv")
