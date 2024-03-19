
##LUAD
library(Seurat)
# Read csv
data1 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219674_LJQ-T.counts.tsv", header = T,row.names = 1,sep = "\t")
data2 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219675_GBG-T.counts.tsv", header = T,row.names= 1,sep = "\t")                                     
data3 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219676_LYB-T1.counts.tsv", header = T,row.names= 1,sep = "\t")
data4 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219677_LYB-T2.counts.tsv", header = T,row.names= 1,sep = "\t")
data5 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219678_CYD-T.counts.tsv", header = T,row.names= 1,sep = "\t")
data6 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219679_CYZ-T.counts.tsv", header = T,row.names= 1,sep = "\t")
data7 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219680_XMS.counts.tsv", header = T,row.names= 1,sep = "\t")
data8 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219681_ZYQ.counts.tsv", header = T,row.names= 1,sep = "\t")
data9 <- read.csv("GSE171145_RAW/GSE171145_RAW/GSM5219682_TGS.counts.tsv", header = T,row.names= 1,sep = "\t")

{
  seurat_obj1 <- CreateSeuratObject(counts = data1,project ="GSM5219674" )
seurat_obj2 <- CreateSeuratObject(counts = data2,project = "GSM5219675")
seurat_obj3 <- CreateSeuratObject(counts = data3,project = "GSM5219676")
seurat_obj4 <- CreateSeuratObject(counts = data4,project = "GSM5219677")
seurat_obj5 <- CreateSeuratObject(counts = data5,project = "GSM5219678")
seurat_obj6 <- CreateSeuratObject(counts = data6,project = "GSM5219679")
seurat_obj7 <- CreateSeuratObject(counts = data7,project = "GSM5219680")
seurat_obj8 <- CreateSeuratObject(counts = data8,project = "GSM5219681")
seurat_obj9 <- CreateSeuratObject(counts = data9,project = "GSM5219682")}

scRNA_LUAD <- merge(seurat_obj1, 
                   y = c(seurat_obj2, seurat_obj3,seurat_obj4, seurat_obj5,seurat_obj6, seurat_obj7,seurat_obj8, seurat_obj9)
)


###BRCA
# read csv
data1 <- read.table("GSE148673_RAW/GSE148673_RAW/GSM4476486_combined_UMIcount_CellTypes_TNBC1.txt", header = T,row.names = 1,sep = "\t",skip = 2)
data2 <- read.table("GSE148673_RAW/GSE148673_RAW/GSM4476487_combined_UMIcount_CellTypes_TNBC2.txt", header = T,row.names = 1,sep = "\t",skip = 2)
data3 <- read.table("GSE148673_RAW/GSE148673_RAW/GSM4476488_combined_UMIcount_CellTypes_TNBC3.txt", header = T,row.names = 1,sep = "\t",skip = 2)
data4 <- read.table("GSE148673_RAW/GSE148673_RAW/GSM4476489_combined_UMIcount_CellTypes_TNBC4.txt", header = T,row.names = 1,sep = "\t",skip = 2)
data5 <- read.table("GSE148673_RAW/GSE148673_RAW/GSM4476490_combined_UMIcount_CellTypes_TNBC5.txt", header = T,row.names = 1,sep = "\t",skip = 2)



{
seurat_obj1 <- CreateSeuratObject(counts = data1,project = "GSM4476486")
seurat_obj2 <- CreateSeuratObject(counts = data2,project = "GSM4476487")
seurat_obj3 <- CreateSeuratObject(counts = data3,project = "GSM4476488")
seurat_obj4 <- CreateSeuratObject(counts = data4,project = "GSM4476489")
seurat_obj5 <- CreateSeuratObject(counts = data5,project = "GSM4476490")
}

scRNA_BRCA <- merge(seurat_obj1, 
                    y = c(seurat_obj2, seurat_obj3,seurat_obj4, seurat_obj5))

#LIHC
#Read10x
dir_name=c('GSM5076749','GSM5076750')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE166635_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x)
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i])
}

scRNA_LIHC <- merge(datalist[[1]],
                    y = datalist[[2]],
                    project =dir_name)
save(scRNA_BRCA,scRNA_LUAD,scRNA_LIHC,file = "scRNA_all.Rdata")

load("scRNA_all.Rdata")

GenesPerUMI
scRNA_BRCA[["percent.mt"]] <- PercentageFeatureSet(scRNA_BRCA, pattern = "^MT-")
scRNA_LIHC[["percent.mt"]] <- PercentageFeatureSet(scRNA_LIHC, pattern = "^MT-")
scRNA_LUAD[["percent.mt"]] <- PercentageFeatureSet(scRNA_LUAD, pattern = "^MT-")

## count log10GenesPerUMI and add it to seurat object metadata
scRNA_BRCA$log10GenesPerUMI <- log10(scRNA_BRCA$nFeature_RNA)/log10(scRNA_BRCA$nCount_RNA)
scRNA_LIHC$log10GenesPerUMI <- log10(scRNA_LIHC$nFeature_RNA)/log10(scRNA_LIHC$nCount_RNA)
scRNA_LUAD$log10GenesPerUMI <- log10(scRNA_LUAD$nFeature_RNA)/log10(scRNA_LUAD$nCount_RNA)

#Calculate the proportion of red blood cells
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA_BRCA@assays$RNA)) 
HB.genes <- rownames(scRNA_BRCA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA_BRCA[["percent.HB"]]<-PercentageFeatureSet(scRNA_BRCA, features=HB.genes) 

HB_m <- match(HB.genes, rownames(scRNA_LIHC@assays$RNA)) 
HB.genes <- rownames(scRNA_LIHC@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA_LIHC[["percent.HB"]]<-PercentageFeatureSet(scRNA_LIHC, features=HB.genes)

HB_m <- match(HB.genes, rownames(scRNA_LUAD@assays$RNA)) 
HB.genes <- rownames(scRNA_LUAD@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA_LUAD[["percent.HB"]]<-PercentageFeatureSet(scRNA_LUAD, features=HB.genes)
#"nFeature_RNA", "nCount_RNA", "percent.mt"
p1<-VlnPlot(scRNA_BRCA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 2, pt.size=0.05)
p2<-VlnPlot(scRNA_LIHC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 2, pt.size=0.05)
p3<-VlnPlot(scRNA_LUAD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 2, pt.size=0.05)
dir.create("OUT")
library(ggplot2)
ggsave("OUT/BRCA_QC_before.pdf",p1,width = 7,height = 7)
ggsave("OUT/LIHC_QC_before.pdf",p2,width = 7,height = 7)
ggsave("OUT/LUAD_QC_before.pdf",p3,width = 7,height = 7)


p1 <- FeatureScatter(scRNA_BRCA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() 
p2 <- FeatureScatter(scRNA_BRCA, feature1 = "nCount_RNA", feature2 = "percent.mt" ) + NoLegend() 
p3 <- FeatureScatter(scRNA_BRCA, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend() 
pearplot <- CombinePlots(plots = list(p1,p2,p3), nrow=1, legend="none") 
ggsave("OUT/BRCA_pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 

p1 <- FeatureScatter(scRNA_LIHC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() 
p2 <- FeatureScatter(scRNA_LIHC, feature1 = "nCount_RNA", feature2 = "percent.mt" ) + NoLegend() 
p3 <- FeatureScatter(scRNA_LIHC, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend() 
pearplot <- CombinePlots(plots = list(p1,p2,p3), nrow=1, legend="none") 
ggsave("OUT/LIHC_pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5)

p1 <- FeatureScatter(scRNA_LUAD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() 
p2 <- FeatureScatter(scRNA_LUAD, feature1 = "nCount_RNA", feature2 = "percent.mt" ) + NoLegend() 
p3 <- FeatureScatter(scRNA_LUAD, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend() 
pearplot <- CombinePlots(plots = list(p1,p2,p3), nrow=1, legend="none") 
ggsave("OUT/LUAD_pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5)

table(grepl("^ERCC-",rownames(scRNA_BRCA)))
table(grepl("^ERCC-",rownames(scRNA_LIHC)))
table(grepl("^ERCC-",rownames(scRNA_LUAD)))
#nUMI > 500
#nGene > 250
#log10GenesPerUMI > 0.8
#mitoRatio < 0.2


scRNA_BRCA <- subset(x = scRNA_BRCA, subset= (nFeature_RNA >= 200) & (nFeature_RNA <= 6000) &
                          (log10GenesPerUMI > 0.80) & (percent.HB <= 1) &
                          (percent.mt < 20))   
scRNA_LIHC <- subset(x =scRNA_LIHC, subset= (nFeature_RNA >= 200) & (nFeature_RNA <= 6000) &
                       (log10GenesPerUMI > 0.80) & (percent.HB <= 1) &
                       (percent.mt < 20))  
scRNA_LUAD <- subset(x = scRNA_LUAD, subset= (nFeature_RNA >= 200) & (nFeature_RNA <= 6000) &
                       (log10GenesPerUMI > 0.80) & (percent.HB <= 1) &
                       (percent.mt < 20))  





#Clustering analysis was performed for each tumor
scRNA<-scRNA_LIHC

scRNA<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
scRNA<-FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000)
top10<-head(VariableFeatures(scRNA),10)
plot1<-VariableFeaturePlot(scRNA)
plot2<-LabelPoints(plot = plot1,points = top10,repel = TRUE)
dev.new()
p<-plot1+plot2
ggsave("OUT/variable_LUAD.pdf", plot = p, width = 12, height = 5)

# Scale data
all.genes<-rownames(scRNA)
scRNA<-ScaleData(scRNA,features = all.genes)
#save(scRNA,file = "scRNAscale.Rdata")
save(scRNA,file = "output/scale_LUAD.Rdata")

scRNA<-RunPCA(scRNA,features = VariableFeatures(object=scRNA),seed.use = 1)
print(scRNA[["pca"]],dims = 1:5,nfeatures = 5)

VizDimLoadings(scRNA,dims = 1:2,reduction = "pca")
DimPlot(scRNA,reduction = "pca")
DimHeatmap(scRNA,dims = 1,cells = 500,balanced = TRUE)

# Determine the dimensions of the dataset
scRNA<-JackStraw(scRNA,num.replicate = 100)
scRNA<-ScoreJackStraw(scRNA,dims = 1:20)

JackStrawPlot(scRNA,dims = 1:15)
ElbowPlot(scRNA,ndims = 50)

# Clustered cells
scRNA<-FindNeighbors(scRNA,dims = 1:30)
scRNA_cluster <- scRNA
library(clustree)
scRNA_cluster <- FindClusters(
  object = scRNA_cluster,
  resolution = c(seq(.2,1.2,.1))
)
p<-clustree(scRNA_cluster@meta.data, prefix = "RNA_snn_res.")

scRNA<-FindClusters(scRNA,resolution = 0.7)
head(Idents(scRNA),5)

# Run nonlinear dimensionality reduction（UMAP/tSNE）
scRNA<-RunUMAP(scRNA,dims = 1:30,seed.use = 1)
DimPlot(scRNA,reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)
library(dplyr)

#load("scRNA_luad.Rdata")

scRNA.markers<-FindAllMarkers(scRNA,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 )
  #dir.create("Files")
LUAD_markers<-scRNA.markers
LUAD_top10<-scRNA.markers %>%
  group_by(cluster) %>%
  slice_max(n=10,order_by = avg_log2FC)
scRNA_LUAD <- scRNA
#save markers
dir.create("output")
write.csv(LUAD_top10,file = "output/LUAD_10.csv")
write.csv(LUAD_markers,file = "output/LUAD_markers.csv")
save(scRNA_LUAD,file = "output/scRNA_LUAD.Rdata")



library(Seurat)
####
#Load breast cancer single-cell data
load("scRNA_BRCA.Rdata")
BRCA.marker<-read.csv("BRCA_markers.csv")
head(BRCA.marker)
dim(BRCA.marker)
#[1] 24344     8
library(tidyverse)
all.markers = BRCA.marker %>% dplyr::select(gene, everything()) %>% subset(BRCA.marker$p_val<0.05 & abs(BRCA.marker$avg_log2FC) > 0.5)
#An adjusted P value < 0.05and | log 2 [fold change (FC)] | > 0.5 
dim(all.markers)
#[1] 11717     8
summary(all.markers)

#Pick the top ten most significant genes for each as a visualization of the heat map
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA_BRCA)) 
top10
length(top10)
length(unique(sort(top10)))

#Heat maps are drawn and visualized on the first ten markers of differential expressions
p <- DoHeatmap(scRNA_BRCA, features = unique(top10))+scale_fill_gradientn(colors = c("#0068D2", "white", "tomato"))
ggsave("OUT/pheatmap_BRCA_marker.pdf",p,width = 15,height = 28)

p<-DotPlot(scRNA_BRCA, features = unique(top10))+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1 ,size = 8))
ggsave("OUT/dotplot_BRCA_marker.pdf",p,width = 28,height = 15)



library(SingleR)
refdata <- get(load("HumanPrimaryCellAtlasData.Rdata"))
assay(refdata)[1:4,1:4]
head(refdata@colData)
head(refdata)
library(celldex)
library(Seurat)
testdata <- GetAssayData(scRNA_BRCA, slot="data")
clusters <- scRNA_BRCA@meta.data$RNA_snn_res.0.7
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    #labels = hpca.se$label.main,#labels = hpca.se$label.fine
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)


table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) 
#Combined with the above results, celltype annotation information was added to the scRNA
scRNA_BRCA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_BRCA@meta.data[which(scRNA_BRCA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 <- DimPlot(scRNA_BRCA, group.by="celltype", label=F )
p1
ggsave("OUT/celltype_BRCA.pdf",p1,width = 10,height = 6)
save(scRNA_BRCA,file = "Files/scRNA_BRCA_celltype.Rdata")


rm(list=ls())
####
#Load lung cancer single-cell data
load("scRNA_LUAD.Rdata")
LUAD.marker<-read.csv("LUAD_markers.csv")
head(LUAD.marker)
dim(LUAD.marker)
#[1] 25913     8
library(tidyverse)
all.markers = LUAD.marker %>% dplyr::select(gene, everything()) %>% subset(LUAD.marker$p_val<0.05 & abs(LUAD.marker$avg_log2FC) > 0.5)
dim(all.markers)
#[1] 13556     8
summary(all.markers)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA_LUAD)) 
top10
length(top10)
length(unique(sort(top10)))

p <- DoHeatmap(scRNA_LUAD, features = unique(top10))+scale_fill_gradientn(colors = c("#0068D2", "white", "tomato"))
ggsave("OUT/pheatmap_LUAD_marker.pdf",p,width = 15,height = 28)
p<-DotPlot(scRNA_LUAD, features = unique(top10))+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1 ,size = 8))
ggsave("OUT/dotplot_LUAD_marker.pdf",p,width = 30,height = 13)




library(SingleR)
refdata <- get(load("HumanPrimaryCellAtlasData.Rdata"))
assay(refdata)[1:4,1:4]
head(refdata@colData)
head(refdata)
library(celldex)
library(Seurat)
testdata <- GetAssayData(scRNA_LUAD, slot="data")
clusters <- scRNA_LUAD@meta.data$RNA_snn_res.0.7
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)


table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) 
scRNA_LUAD@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_LUAD@meta.data[which(scRNA_LUAD@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 <- DimPlot(scRNA_LUAD, group.by="celltype", label=F )
p1
ggsave("OUT/celltype_LUAD2.pdf",p1,width = 10,height = 6)
save(scRNA_LUAD,file = "Files/scRNA_LUAD_celltype.Rdata")



rm(list=ls())

load("scRNA_LIHC.Rdata")
LIHC.marker<-read.csv("LIHC_markers.csv")
head(LIHC.marker)
dim(LIHC.marker)
#17092     8
library(tidyverse)
all.markers = LIHC.marker %>% dplyr::select(gene, everything()) %>% subset(LIHC.marker$p_val<0.05 & abs(LIHC.marker$avg_log2FC) > 0.5)
dim(all.markers)
summary(all.markers)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA_LIHC)) 
top10
length(top10)
length(unique(sort(top10)))

p <- DoHeatmap(scRNA_LIHC, features = unique(top10))+scale_fill_gradientn(colors = c("#0068D2", "white", "tomato"))
ggsave("OUT/pheatmap_LIHC_marker.pdf",p,width = 15,height = 28)
p<-DotPlot(scRNA_LIHC, features = unique(top10))+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1 ,size = 8))
ggsave("OUT/dotplot_LIHC_marker.pdf",p,width = 30,height = 13)




library(SingleR)
assay(refdata)[1:4,1:4]
head(refdata@colData)
head(refdata)
library(celldex)
library(Seurat)
testdata <- GetAssayData(scRNA_LIHC, slot="data")
clusters <- scRNA_LIHC@meta.data$RNA_snn_res.0.5
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)


table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) 
scRNA_LIHC@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_LIHC@meta.data[which(scRNA_LIHC@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 <- DimPlot(scRNA_LIHC, group.by="celltype", label=F )
p1

FeaturePlot(scRNA,features = 'SPP1')
VlnPlot(scRNA, 
        features = 'SPP1',group.by = 'celltype',
        pt.size = 0) 
save(scRNA,file = "scRNA_LIHC.Rdata")



rm(list=ls())
library(Seurat)
library(gplots)
library(ggplot2)
library(dplyr)
##Data from all macrophages were extracted

########
#BRCA
load("Files/scRNA_BRCA_celltype.Rdata")
Macro_BRCA = scRNA_BRCA[,scRNA_BRCA@meta.data$celltype %in% "Macrophage"]
save(Macro_BRCA,file = "Files/Macro.Rdata")

########
#LIHC
load("Files/scRNA_LIHC_celltype.Rdata")
Macro_LIHC = scRNA_LIHC[,scRNA_LIHC@meta.data$celltype %in% "Macrophage"]
load("Files/Macro.Rdata")
save(Macro_BRCA,Macro_LIHC,file = "Files/Macro.Rdata")


########
#LUAD
load("Files/scRNA_LUAD_celltype.Rdata")
Macro_LUAD = scRNA_LUAD[,scRNA_LUAD@meta.data$celltype %in% "Macrophage"]
load("Files/Macro.Rdata")
save(Macro_BRCA,Macro_LIHC,Macro_LUAD,file = "Files/Macro.Rdata")

#########
rm(list=ls())
load("Files/Macro.Rdata")
#####LIHC
p1<-DimPlot(Macro_LIHC, group.by="seurat_clusters", label=F )
ggsave("OUT/macro_LIHC_dimplot.pdf",p1,width = 5,height = 4)

LIHC_markers<-FindAllMarkers(Macro_LIHC,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 )
save(LIHC_markers, file = "Files/Macro_marker.Rdata")
markers = LIHC_markers %>% select(gene, everything()) %>% subset(LIHC_markers$p_val<0.05 & abs(LIHC_markers$avg_log2FC) > 0.5)
p2<-DoHeatmap(Macro_LIHC, label = F,size = 2,features = unique(markers$gene))+scale_fill_gradientn(colors = c("#0068D2", "white", "tomato"))
ggsave("OUT/macro_LIHC_heatmap.pdf",p2,width = 6,height = 10)

top10 = markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(Macro_LIHC)) 
p<-DotPlot(Macro_LIHC, features = unique(top10))+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1 ,size = 8))
ggsave("OUT/macro_LIHC_dotplot.pdf",p,width = 10,height = 5)



#####BRCA
p1<-DimPlot(Macro_BRCA, group.by="seurat_clusters", label=F )
ggsave("OUT/macro_BRCA_dimplot.pdf",p1,width = 5,height = 4)

BRCA_markers<-FindAllMarkers(Macro_BRCA,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 )
save(BRCA_markers,LIHC_markers, file = "Files/Macro_marker.Rdata")
markers = BRCA_markers %>%  dplyr ::select(gene, everything()) %>% subset(BRCA_markers$p_val<0.05 & abs(BRCA_markers$avg_log2FC) > 0.5)
p2<-DoHeatmap(Macro_BRCA, label = F,size = 2,features = unique(markers$gene))+scale_fill_gradientn(colors = c("#0068D2", "white", "tomato"))
ggsave("OUT/macro_BRCA_heatmap.pdf",p2,width = 6,height = 10)
top10 = markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(Macro_BRCA)) 
p<-DotPlot(Macro_BRCA, features = unique(top10))+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1 ,size = 8))

#####LUAD
p1<-DimPlot(Macro_LUAD, group.by="seurat_clusters", label=F )
ggsave("OUT/macro_LUAD_dimplot.pdf",p1,width = 5,height = 4)

LUAD_markers<-FindAllMarkers(Macro_LUAD,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 )
save(LUAD_markers,LIHC_markers,BRCA_markers, file = "Files/Macro_marker.Rdata")
markers = LUAD_markers %>%  dplyr ::select(gene, everything()) %>% subset(LUAD_markers$p_val<0.05 & abs(LUAD_markers$avg_log2FC) > 0.5)
p2<-DoHeatmap(Macro_LUAD, label = F,size = 2,features = unique(markers$gene))+scale_fill_gradientn(colors = c("#0068D2", "white", "tomato"))
ggsave("OUT/macro_LUAD_heatmap.pdf",p2,width = 6,height = 10)

top10 = markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(Macro_LUAD)) 
p<-DotPlot(Macro_LUAD, features = unique(top10))+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1 ,size = 8))
ggsave("OUT/macro_LUAD_dotplot.pdf",p,width = 10,height = 5)


######KEGG
rm(list = ls())
load("Files/Macro_marker.Rdata")
library(clusterProfiler)
library(org.Hs.eg.db)
BRCA_markers<-BRCA_markers[which(BRCA_markers$p_val<0.05),]
ids=bitr(BRCA_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## SYMBOL to ENTREZID
BRCA_markers=merge(BRCA_markers,ids,by.x='gene',by.y='SYMBOL')

gcSample=split(BRCA_markers$ENTREZID, BRCA_markers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p<-p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("OUT/BRCA_KEGG.pdf",p,width = 5,height = 6)

LIHC_markers<-LIHC_markers[which(LIHC_markers$p_val<0.05),]
ids=bitr(LIHC_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## SYMBOL to ENTREZID
LIHC_markers=merge(LIHC_markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(LIHC_markers$ENTREZID, LIHC_markers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p <-p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("OUT/LIHC_KEGG.pdf",p,width = 5,height = 6)
LUAD_markers<-LUAD_markers[which(LUAD_markers$p_val<0.05),]
ids=bitr(LUAD_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ##  SYMBOL to ENTREZID
LUAD_markers=merge(LUAD_markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(LUAD_markers$ENTREZID, LUAD_markers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p<-p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("OUT/LUAD_KEGG.pdf",p,width = 5,height = 6)
