##   Correlation between SPP1 and immune cells in different tumors

exprs<-list()

library(TCGAbiolinks)
library(tidyr)
projects <- getGDCprojects()
library(dplyr)
projects <- projects %>% 
  as.data.frame() %>% 
  dplyr::select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))
a<-c("LIHC","CESC","LUAD","COAD","BRCA","ESCA","KIRP","STAD","PRAD","UCEC","HNSC","KIRC","LUSC","THCA","READ","BLCA","KICH","CHOL")
b<-match(a,projects$tumor)
projects<-projects[b,]




for (i in 1:18) {
  filename1 = paste0("KM/result/",projects$project_id[i],"-tpm.txt")
  mRNA<-read.table(filename1,header = T,check.names = F)
  mRNA<-mRNA[,-c(1,3)]
  data <- aggregate(x = mRNA[,-1],
                    by = list(mRNA$gene_name),
                    FUN = mean)
  gene_name<-data[,1]
  logTPM <- log2(data[,-1]+1)
  #boxplot(logTPM[,1:10],las=2)
  #head(logTPM[,1:3])
  group_list=ifelse(as.numeric(substr(colnames(logTPM),14,15)) < 10,'tumor','normal')
  group_list<-data.frame(group_list)
  group_list$ID<-colnames(logTPM)
  group_list<-group_list[which(group_list$group_list=="tumor"),]
  logTPM<-logTPM[,group_list$ID]
  logTPM$gene_name<-gene_name
  exprs[[i]]<-logTPM
}

save(exprs,file = "exprs.Rdata")

##########
##ssGSEA
load("exprs.Rdata")
{
  library(genefilter)
  library(GSVA)
  library(Biobase)
  library(stringr)
  library(GSEABase)
  library(pheatmap)
}

for (i in 1:18) {
#Read a single tumor matrix
tumor<-as.data.frame(exprs[[i]])
tumor[1:7,1:4]
rownames(tumor)<-tumor$gene_name
tumor<-dplyr::select(tumor,-gene_name)
gene_set<-read.csv("mmc.csv")[, 1:2] 
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#SPP1
SPP1<-data.frame(tumor[which(rownames(tumor)=="SPP1"),])
SPP1<-data.frame(t(SPP1))
SPP1$ID<-rownames(SPP1)
SPP1<-data.frame(SPP1[order(SPP1$SPP1,decreasing = F),])#Expression data is sorted
SPP1$ID<-gsub("[.]","-",SPP1$ID)#If the sample name changes during the matrix conversion process, the name is restored
#计算ssGSEA
#dir.create("ssGSEA")
gsva_matrix<- gsva(as.matrix(tumor),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
filenames<-paste0("ssGSEA/gsva_matrix_",projects$project_id[i],".Rdata")
save(gsva_matrix,file = filenames)
#dir.create("SPP1")
filenames<-paste0("SPP1/SPP1_",projects$project_id[i],".Rdata")

save(SPP1,file = filenames)
}





library(ComplexHeatmap)
library(circlize)
p<-list()
for (i in 1:18) {
filenames<-paste0("ssGSEA/gsva_matrix_",projects$project_id[i],".Rdata")  
load(filenames)
filenames<-paste0("SPP1/SPP1_",projects$project_id[i],".Rdata")
load(filenames)
gsva_matrix<-gsva_matrix[,SPP1$ID]#Samples were arranged in the order in which SPP1 was expressed
gsva_matrix1<- t(scale(t(gsva_matrix))) 
#Define standardized calculations and standardize data
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1)
nor_gsva_matrix1[1:5,1:5]

#Define pro-tumor and anti-tumor genes
anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
anti<- gsub('^ ','',rownames(nor_gsva_matrix1))%in%anti_tumor
pro<- gsub('^ ','',rownames(nor_gsva_matrix1))%in%pro_tumor
non <- !(anti|pro)
nor_gsva_matrix1<- rbind(nor_gsva_matrix1[anti,],nor_gsva_matrix1[pro,],nor_gsva_matrix1[non,])


spp1_col = colorRamp2(c(0,max(SPP1$SPP1)/2,max(SPP1$SPP1)), c("#B0A2A2","#7D6E8F" , "#5B4C6E"))
top_annotation = HeatmapAnnotation(SPP1 = SPP1$SPP1,
                                   col =(list(SPP1 = spp1_col)),
                                   show_legend = T)
p[[i]]<-ComplexHeatmap::pheatmap(nor_gsva_matrix1,
                         show_colnames = F,
                         cluster_rows =F,
                         cluster_cols = F,
                         color =colorRamp2(c(min(nor_gsva_matrix1), max(nor_gsva_matrix1)/2, max(nor_gsva_matrix1)), c("blue", "white", "red")) ,
                         top_annotation = top_annotation, 
                         legend = T,
                         border_color = NA,
                         cellheight = 12,
                         fontsize_row = 10,
                         main = projects$project_id[i])


}


file<-paste0("OUT/ssGSEA/",projects$tumor[1],".png")
Cairo::CairoPNG( 
  filename = file, 
  width = 12,           
  height = 12,          
  units = "in",       
  dpi = 600)          

p[[i]]
dev.off()





#correlation
library("Hmisc")
library(corrplot)
library(Cairo)
for (i in 1:18) {
  
  filenames<-paste0("ssGSEA/gsva_matrix_",projects$project_id[i],".Rdata")  
  load(filenames)
  filenames<-paste0("SPP1/SPP1_",projects$project_id[i],".Rdata")
  load(filenames)
ssgsva<-data.frame(t(gsva_matrix))
ssgsva$ID<-rownames(ssgsva)
ssgsva$ID[1:4]
matrix<-merge(SPP1,ssgsva,by = "ID")

colnames(matrix)[1]<-"SPP1"
matrix[1:5,1:5]

res <- cor(matrix[,-1])
A<-round(res, 2)
res2 <- rcorr(as.matrix(matrix[,-1]))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )}
}


library("Hmisc")
library(corrplot)
library(Cairo)
SPP1cor <- list()
P <- list()
for (i in 1:18) {
  
  filenames<-paste0("ssGSEA/gsva_matrix_",projects$project_id[i],".Rdata")  
  load(filenames)
  filenames<-paste0("SPP1/SPP1_",projects$project_id[i],".Rdata")
  load(filenames)
  ssgsva<-data.frame(t(gsva_matrix))
  ssgsva$ID<-rownames(ssgsva)
  ssgsva$ID[1:4]
  matrix<-merge(SPP1,ssgsva,by = "ID")
  
  colnames(matrix)[1]<-"SPP1"
  matrix[1:5,1:5]
  
  res <- cor(matrix[,-1])
  A<-round(res, 2)
  res2 <- rcorr(as.matrix(matrix[,-1]))
  data <- as.data.frame(res2$r)
  SPP1cor[[i]] <- data$SPP1
  pvalue <- as.data.frame(res2$P)
  P[[i]] <- pvalue$SPP1

  }
data <- do.call(rbind,SPP1cor)
colnames(data) <- colnames(res)
data[,1] <- projects[,2]
rownames(data) <- data[,1]
data <- data[,-1]
data1<- as.data.frame(lapply(data,function(x) as.numeric(as.character(x))))
rownames(data1) <- rownames(data)
pvalue <- do.call(rbind,P)
colnames(pvalue) <- colnames(res)
pvalue[,1] <- projects[,2]
library(pheatmap)
bk <- c(seq(-0.8,-0.1,by=0.01),seq(0,0.8,by=0.01))
pheatmap(data1,
         cluster_rows = F,
         breaks = bk,
         color = colorRampPalette(c("blue", "white", "red"))(length(bk)))

figP <- ifelse(pvalue>0.05,"NS",ifelse(pvalue>0.01,"*",ifelse(pvalue>0.001,"**",ifelse(pvalue>0.001,"***","****"))))
figP <- figP[,-1]
pheatmap(data1,
         annotation_row = NULL,
         annotation_names_row = T,
         labels_row = rownames(data1),
         legend = T,
         show_rownames = T, 
         angle_col = 90,
         border_color = "#EEEEEE",
         na_col = "white",
         treeheight_row = 0, 
         treeheight_col = 0,
         legend_labels = "sig*log(p-value)",
         color = colorRampPalette(c("blue", "white", "red"))(length(bk)),
         display_numbers = figP,
         breaks = bk,
         filename = 'corpancancer.png')



## Single-cell data reads and expression of SPP1 in epithelial cells and macrophages  ##

library(Seurat)
dir_name=c('GSM7290763','GSM7290769','GSM7290772','GSM7290773','GSM7290774','GSM7290777')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE231559_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x)
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i])
}
scRNA_COAD <- merge(datalist[[1]],
                    y = c(datalist[[2]],datalist[[3]],datalist[[4]],datalist[[5]],datalist[[6]]),
                    project =dir_name)

scRNA_COAD[["percent.mt"]] <- PercentageFeatureSet(scRNA_COAD, pattern = "^MT-")
scRNA_COAD$log10GenesPerUMI <- log10(scRNA_COAD$nFeature_RNA)/log10(scRNA_COAD$nCount_RNA)

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA_COAD@assays$RNA)) 
HB.genes <- rownames(scRNA_COAD@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA_COAD[["percent.HB"]]<-PercentageFeatureSet(scRNA_COAD, features=HB.genes) 
VlnPlot(scRNA_COAD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 2, pt.size=0.05)
FeatureScatter(scRNA_COAD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() 
FeatureScatter(scRNA_COAD, feature1 = "nCount_RNA", feature2 = "percent.mt" ) + NoLegend() 
FeatureScatter(scRNA_COAD, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend() 
table(grepl("^ERCC-",rownames(scRNA_COAD)))
scRNA_COAD <- subset(x = scRNA_COAD, subset= (nFeature_RNA >= 200) & (nFeature_RNA <= 8000) &
                       (log10GenesPerUMI > 0.80) & (percent.HB <= 1) &
                       (percent.mt < 20))   



scRNA<-scRNA_COAD
scRNA<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
scRNA<-FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000)


all.genes<-rownames(scRNA)
scRNA<-ScaleData(scRNA,features = all.genes)


scRNA<-RunPCA(scRNA,features = VariableFeatures(object=scRNA),seed.use = 1)
print(scRNA[["pca"]],dims = 1:5,nfeatures = 5)

VizDimLoadings(scRNA,dims = 1:2,reduction = "pca")
DimPlot(scRNA,reduction = "pca")
DimHeatmap(scRNA,dims = 1,cells = 500,balanced = TRUE)
scRNA<-JackStraw(scRNA,num.replicate = 100)
scRNA<-ScoreJackStraw(scRNA,dims = 1:20)

JackStrawPlot(scRNA,dims = 1:15)
ElbowPlot(scRNA,ndims = 50)

scRNA<-FindNeighbors(scRNA,dims = 1:20)
scRNA_cluster <- scRNA
library(clustree)
scRNA_cluster <- FindClusters(
  object = scRNA_cluster,
  resolution = c(seq(.2,1.2,.1))
)
p<-clustree(scRNA_cluster@meta.data, prefix = "RNA_snn_res.")

scRNA<-FindClusters(scRNA,resolution = 0.7)
head(Idents(scRNA),5)

#（UMAP/tSNE）
scRNA<-RunUMAP(scRNA,dims = 1:20,seed.use = 1)
DimPlot(scRNA,reduction = "umap")

#  SingleR 
library(SingleR)
load("HumanPrimaryCellAtlasData.Rdata")
assay(ref)[1:4,1:4]
head(ref@colData)
library(celldex)
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$RNA_snn_res.0.7
cellpred <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)


table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) 
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 <- DimPlot(scRNA, group.by="celltype", label=F ,cols = c(rep('grey',9),'#87CEFA',rep('grey',7),'#9F79EE','grey','pink',rep('grey',13)))
p1

FeaturePlot(scRNA,features = 'SPP1')
VlnPlot(scRNA, 
        features = 'SPP1',group.by = 'celltype',
        pt.size = 0) 
save(scRNA,file = "scRNA_COAD.Rdata")

######
#THCA

library(Seurat)
library(data.table)

data1 <- fread("GSE232237_RAW/GSM7324272_PT9.count.tsv.gz", header = T,sep = "\t")
data2 <- fread("GSE232237_RAW/GSM7324271_PT8.count.tsv.gz", header = T,sep = "\t")
data3 <- fread("GSE232237_RAW/GSM7324270_PT7.count.tsv.gz", header = T,sep = "\t")
data4 <- fread("GSE232237_RAW/GSM7324269_PT5.count.tsv.gz", header = T,sep = "\t")
data5 <- fread("GSE232237_RAW/GSM7324268_PT3.count.tsv.gz", header = T,sep = "\t")
data6 <- fread("GSE232237_RAW/GSM7324267_PT12.count.tsv.gz", header = T,sep = "\t")
data7 <- fread("GSE232237_RAW/GSM7324266_PT10.count.tsv.gz", header = T,sep = "\t")

library(tibble)
data1 <- data1 %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'gene_name')
data2 <- data2 %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'gene_name')
data3 <- data3 %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'gene_name')
data4 <- data4 %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'gene_name')
data5 <- data5 %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'gene_name')
data6 <- data6 %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'gene_name')
data7 <- data7 %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = 'gene_name')


{
  seurat_obj1 <- CreateSeuratObject(counts = data1,project ="GSM7324272" )
  seurat_obj2 <- CreateSeuratObject(counts = data2,project = "GSM7324271")
  seurat_obj3 <- CreateSeuratObject(counts = data3,project = "GSM7324270")
  seurat_obj4 <- CreateSeuratObject(counts = data4,project = "GSM7324269")
  seurat_obj5 <- CreateSeuratObject(counts = data5,project = "GSM7324268")
  seurat_obj6 <- CreateSeuratObject(counts = data6,project = "GSM7324267")
  seurat_obj7 <- CreateSeuratObject(counts = data7,project = "GSM7324266")
  }

scRNA_THCA <- merge(seurat_obj1, 
                    y = c(seurat_obj2, seurat_obj3,seurat_obj4, seurat_obj5,seurat_obj6, seurat_obj7)
)



scRNA_THCA[["percent.mt"]] <- PercentageFeatureSet(scRNA_THCA, pattern = "^MT-")
scRNA_THCA$log10GenesPerUMI <- log10(scRNA_THCA$nFeature_RNA)/log10(scRNA_THCA$nCount_RNA)

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA_THCA@assays$RNA)) 
HB.genes <- rownames(scRNA_THCA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA_THCA[["percent.HB"]]<-PercentageFeatureSet(scRNA_THCA, features=HB.genes) 
VlnPlot(scRNA_THCA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 2, pt.size=0.05)
FeatureScatter(scRNA_THCA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() 
FeatureScatter(scRNA_THCA, feature1 = "nCount_RNA", feature2 = "percent.mt" ) + NoLegend() 
FeatureScatter(scRNA_THCA, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend() 
table(grepl("^ERCC-",rownames(scRNA_THCA)))
scRNA_THCA <- subset(x = scRNA_THCA, subset= (nFeature_RNA >= 200) & (nFeature_RNA <= 8000) &
                       (log10GenesPerUMI > 0.80) & (percent.HB <= 1) &
                       (percent.mt < 20))   

scRNA<-scRNA_THCA
scRNA<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
scRNA<-FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000)

all.genes<-rownames(scRNA)
scRNA<-ScaleData(scRNA,features = all.genes)

scRNA<-RunPCA(scRNA,features = VariableFeatures(object=scRNA),seed.use = 1)
print(scRNA[["pca"]],dims = 1:5,nfeatures = 5)

VizDimLoadings(scRNA,dims = 1:2,reduction = "pca")
DimPlot(scRNA,reduction = "pca")
DimHeatmap(scRNA,dims = 1,cells = 500,balanced = TRUE)

scRNA<-JackStraw(scRNA,num.replicate = 100)
scRNA<-ScoreJackStraw(scRNA,dims = 1:20)

JackStrawPlot(scRNA,dims = 1:15)
ElbowPlot(scRNA,ndims = 50)

# Clustered cells
# A KNN graph is established, and the edge weights between any two elements are refined based on the shared overlap in their local domains
scRNA<-FindNeighbors(scRNA,dims = 1:30)
#save(scRNA,file = "scRNA_neigbbors")
scRNA_cluster <- scRNA
library(clustree)
scRNA_cluster <- FindClusters(
  object = scRNA_cluster,
  resolution = c(seq(.2,1.2,.1))
)
p<-clustree(scRNA_cluster@meta.data, prefix = "RNA_snn_res.")

scRNA<-FindClusters(scRNA,resolution = 0.9)
head(Idents(scRNA),5)

#（UMAP/tSNE）
scRNA<-RunUMAP(scRNA,dims = 1:30,seed.use = 1)
DimPlot(scRNA,reduction = "umap")

#  SingleR 
library(SingleR)
load("HumanPrimaryCellAtlasData.Rdata")
assay(ref)[1:4,1:4]
head(ref@colData)
library(celldex)
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$RNA_snn_res.0.7
cellpred <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)


table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) 
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 <- DimPlot(scRNA, group.by="celltype", label=F )
p1

FeaturePlot(scRNA,features = 'SPP1')
VlnPlot(scRNA, 
        features = 'SPP1',group.by = 'celltype',
        pt.size = 0) 
save(scRNA,file = "scRNA_COAD.Rdata")



################
#PRAD
library(data.table)
file1 <- fread('GSM5353248_PA_PR5269_4_S28_L002_dge.txt.gz')
file2 <- fread('GSM5353247_PA_PR5269_3_S27_L002_dge.txt.gz')
file3 <- fread('GSM5353246_PA_PR5269_2_S26_L002_dge.txt.gz')
file4 <- fread('GSM5353245_PA_PR5269_1_S25_L002_dge.txt.gz')
file5 <- fread('GSM5353229_PA_PR5199-640K_Pool_2_S56_L002_dge.txt.gz')
file6 <- fread('GSM5353228_PA_PR5199-640K_Pool_1_3_S108_L004_dge.txt.gz')
file7 <- fread('GSM5353227_PA_PR5199-193K_Pool_1_2_3_S55_L002_dge.txt.gz')
file8 <- fread('GSM5353226_PA_PR5196-2_Pool_1_2_3_S54_L002_dge.txt.gz')
file9 <- fread('GSM5353225_PA_PR5196-1_Pool_1_2_3_S53_L002_dge.txt.gz')
file10 <- fread('GSM5353224_PA_PR5186_Pool_1_2_3_S27_L001_dge.txt.gz')

library(tibble)
file1<- file1 %>%  tibble::column_to_rownames(var = 'GENE')
file2<- file2 %>%  tibble::column_to_rownames(var = 'GENE')
file3<- file3 %>%  tibble::column_to_rownames(var = 'GENE')
file4<- file4 %>%  tibble::column_to_rownames(var = 'GENE')
file5<- file5 %>%  tibble::column_to_rownames(var = 'GENE')
file6<- file6 %>%  tibble::column_to_rownames(var = 'GENE')
file7<- file7 %>%  tibble::column_to_rownames(var = 'GENE')
file8<- file8 %>%  tibble::column_to_rownames(var = 'GENE')
file9<- file9 %>%  tibble::column_to_rownames(var = 'GENE')
file10<- file10 %>%  tibble::column_to_rownames(var = 'GENE')

library(Seurat)
seurat_obj1 <- CreateSeuratObject(counts = file1,project ="GSM5353248" )
seurat_obj2 <- CreateSeuratObject(counts = file2,project ="GSM5353247" )
seurat_obj3 <- CreateSeuratObject(counts = file3,project ="GSM5353246" )
seurat_obj4 <- CreateSeuratObject(counts = file4,project ="GSM5353245" )
seurat_obj5 <- CreateSeuratObject(counts = file5,project ="GSM5353229" )
seurat_obj6 <- CreateSeuratObject(counts = file6,project ="GSM5353228" )
seurat_obj7 <- CreateSeuratObject(counts = file7,project ="GSM5353227" )
seurat_obj8 <- CreateSeuratObject(counts = file8,project ="GSM5353226" )
seurat_obj9 <- CreateSeuratObject(counts = file9,project ="GSM5353225" )
seurat_obj10 <- CreateSeuratObject(counts = file10,project ="GSM5353224" )

scRNA_PRAD <- merge(seurat_obj1, 
                    y = c(seurat_obj2, seurat_obj3,seurat_obj4, seurat_obj5,seurat_obj6, seurat_obj7,seurat_obj8,seurat_obj9,seurat_obj10))
                    

scRNA_PRAD[["percent.mt"]] <- PercentageFeatureSet(scRNA_PRAD, pattern = "^MT-")
scRNA_PRAD$log10GenesPerUMI <- log10(scRNA_PRAD$nFeature_RNA)/log10(scRNA_PRAD$nCount_RNA)

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA_PRAD@assays$RNA)) 
HB.genes <- rownames(scRNA_PRAD@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA_PRAD[["percent.HB"]]<-PercentageFeatureSet(scRNA_PRAD, features=HB.genes) 
VlnPlot(scRNA_PRAD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 2, pt.size=0.05,raster=FALSE)
FeatureScatter(scRNA_PRAD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() 
FeatureScatter(scRNA_PRAD, feature1 = "nCount_RNA", feature2 = "percent.mt" ) + NoLegend() 
FeatureScatter(scRNA_PRAD, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend() 
table(grepl("^ERCC-",rownames(scRNA_PRAD)))
scRNA_PRAD <- subset(x = scRNA_PRAD, subset= (nFeature_RNA >= 200) & (nFeature_RNA <= 8000) &
                       (log10GenesPerUMI > 0.80) & (percent.HB <= 1) &
                       (percent.mt < 20))   

scRNA<-scRNA_PRAD
scRNA<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
scRNA<-FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000)

library(Seurat)
all.genes<-rownames(scRNA)
scRNA<-ScaleData(scRNA,features = all.genes)

scRNA<-RunPCA(scRNA,features = VariableFeatures(object=scRNA),seed.use = 1)
print(scRNA[["pca"]],dims = 1:5,nfeatures = 5)

VizDimLoadings(scRNA,dims = 1:2,reduction = "pca")
DimPlot(scRNA,reduction = "pca")
DimHeatmap(scRNA,dims = 1,cells = 500,balanced = TRUE)

scRNA<-JackStraw(scRNA,num.replicate = 100)
scRNA<-ScoreJackStraw(scRNA,dims = 1:20)

JackStrawPlot(scRNA,dims = 1:15)
ElbowPlot(scRNA,ndims = 50)


scRNA<-FindNeighbors(scRNA,dims = 1:25)
scRNA_cluster <- scRNA
library(clustree)
scRNA_cluster <- FindClusters(
  object = scRNA_cluster,
  resolution = c(seq(.2,1.2,.1))
)
p<-clustree(scRNA_cluster@meta.data, prefix = "RNA_snn_res.")

scRNA<-FindClusters(scRNA,resolution = 0.7)
head(Idents(scRNA),5)

# （UMAP/tSNE）
scRNA<-RunUMAP(scRNA,dims = 1:25,seed.use = 1)
DimPlot(scRNA,reduction = "umap")

scRNA<-RunTSNE(scRNA,dims = 1:25,seed.use = 1)
DimPlot(scRNA,reduction = "tsne")

library(SingleR)
load("readin/HumanPrimaryCellAtlasData.Rdata")
assay(ref)[1:4,1:4]
head(ref@colData)
library(celldex)
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$RNA_snn_res.0.7
cellpred <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                    method = clusters,clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)


table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) 
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(rownames(scRNA@meta.data) == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}

p1 <- DimPlot(scRNA, group.by="celltype", label=F ,cols = c(rep('grey',7),'#87CEFA',rep('grey',7),'#9F79EE','grey','pink',rep('grey',11)))
p1

FeaturePlot(scRNA,features = 'SPP1')
VlnPlot(scRNA, 
        features = 'SPP1',group.by = 'celltype',
        pt.size = 0) 
save(scRNA,file = "scRNA_PRAD.Rdata")


