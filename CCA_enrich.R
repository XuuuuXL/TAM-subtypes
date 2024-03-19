library(Seurat)

rm(list=ls())
###Start consolidating data
load("E:\\yjs\\bioinformation\\cancer\\cellmarker\\Files\\Macro.Rdata")
#Merging seurat objects using the merge function
scRNAlist <- list()
#The following code will create a seurat object for each sample and store it in the scRNAlist
  scRNAlist[[1]] <- Macro_BRCA
  scRNAlist[[2]] <- Macro_LIHC
  scRNAlist[[3]] <- Macro_LUAD

  names(scRNAlist)[1]<-"BRCA"
  names(scRNAlist)[2]<-"LIHC"
  names(scRNAlist)[3]<-"LUAD"

for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)

scRNA <- IntegrateData(anchorset = scRNA.anchors)
DefaultAssay(scRNA) <- "integrated"
dim(scRNA)
#[1] 2000 7801

#dir.create("Files")
save(scRNA,file = "Files/scRNA_acc.Rdata")
library(ggplot2)
#Dimensionality reduction clustering

scRNA1 <- ScaleData(scRNA, features = VariableFeatures(scRNA))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1),seed.use = 1)
plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA1, ndims=50, reduction="pca") 
plotc <- plot1+plot2
#dir.create("OUT")
ggsave("OUT/pca.pdf", plot = plotc, width = 8, height = 4)
#principal component selection
pc.num=1:30

##Cell clustering
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num)
scRNA_cluster <- scRNA1
library(clustree)
scRNA_cluster <- FindClusters(
  object = scRNA_cluster,
  resolution = c(seq(.2,1.2,.1))
)
p<-clustree(scRNA_cluster@meta.data, prefix = "integrated_snn_res.")

scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
head(Idents(scRNA1),5)

table(scRNA1@meta.data$seurat_clusters)
#  0    1    2    3    4    5    6    7    8    9   10 
# 2738 1300 1299  466  446  425  408  252  208  182   77 
metadata <- scRNA1@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
#write.csv(cell_cluster,'cluster1/cell_cluster.csv',row.names = F)

##Nonlinear dimensionality reduction
#tSNE
scRNA1 = RunTSNE(scRNA1, dims = pc.num,seed.use = 1)
embed_tsne <- Embeddings(scRNA1, 'tsne')  
#write.csv(embed_tsne,'cluster1/embed_tsne.csv')
#group_by_cluster
plot1 = DimPlot(scRNA1, reduction = "tsne", label=T) 
ggsave("OUT/tSNE.pdf", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(scRNA1, reduction = "tsne", group.by='orig.ident') 
ggsave("OUT/tSNE_sample.pdf", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("OUT/tSNE_cluster_sample.pdf", plot = plotc, width = 10, height = 5)

#UMAP
scRNA1 <- RunUMAP(scRNA1, dims = pc.num,seed.use = 1)
embed_umap <- Embeddings(scRNA1, 'umap')   #Extract the coordinates of the UMAP diagram
#write.csv(embed_umap,'cluster1/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
ggsave("OUT/UMAP_cluster.pdf", plot = plot3, width = 8, height = 7)
ggsave("E:\\yjs\\bioinformation\\cancer\\Plots\\p1\\UMAP_sample.png",plot = plot3,width = 8, height = 7)

#group_by_sample
plot4 = DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')
ggsave("OUT/UMAP_sample.pdf", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("OUT/UMAP_cluster_sample.pdf", plot = plotc, width = 10, height = 5)

#Merge tSNE with UMAP
library(patchwork)
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("OUT/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)

save(scRNA1,file = "Files/scRNA_tsne_umap.Rdata")


##==Data quality control==#
scRNA <- scRNA1  #Subsequent analyses were performed using the consolidated data
##meta.data add information
proj_name <- data.frame(proj_name=rep("demo2",ncol(scRNA)))
rownames(proj_name) <- row.names(scRNA@meta.data)
scRNA <- AddMetaData(scRNA, proj_name)

##Switching data sets
DefaultAssay(scRNA) <- "RNA"

##Calculate mitochondrial and red blood cell gene ratios
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#Calculate red blood cell ratio
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))

##Draw a violin diagram
#All samples are used for violin figures group.by="proj_name"，One violin per sample group.by="orig.ident"
violin <-VlnPlot(scRNA, group.by = "proj_name",  
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.01, #You don't need to display the dots, you can set them pt.size = 0
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("OUT/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("OUT/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 

minGene=500
maxGene=5000
pctMT=20
##Data quality control
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))
violin <-VlnPlot(scRNA, group.by = "proj_name",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("OUT/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 

#Transform the dataset
DefaultAssay(scRNA) <- "integrated"
library(dplyr)
scRNA.markers<-FindAllMarkers(scRNA,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 )
# %>% which is a pipeline operation, which means passing the object to the left of %>% to the function on the right, as a setting for the first option (or the setting for the only option left)
all.markers = scRNA.markers %>% subset(scRNA.markers$p_val<0.05 & abs(scRNA.markers$avg_log2FC) > 0.5)
dim(all.markers)
summary(all.markers)

#Pick the top ten most significant genes for each as a visualization of the heat map
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(top10[,c(1,2,5,6,7)],file = "top10.csv")
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA)) 
top10
length(top10)
length(unique(sort(top10)))#There are duplications, i.e., there are genes that are significant in multiple classes

#Heat maps are drawn and visualized on the first ten markers of differential expressions
p <- DoHeatmap(scRNA, features = unique(top10))+scale_fill_gradientn(colors = c("blue", "white", "red"))
ggsave("OUT/pheatmap_marker.pdf",p,width = 18,height = 20)
ggsave("E:\\yjs\\bioinformation\\cancer\\Plots\\p1\\pheatmap_marker.png",p,width = 18,height = 20)
p<-DotPlot(scRNA, features = unique(top10))+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1 ,size = 8))
ggsave("OUT/dotplot_marker.pdf",p,width = 20,height = 7)
ggsave("E:\\yjs\\bioinformation\\cancer\\Plots\\p1\\dotplot_marker.png",p,width = 15,height = 5)
save(scRNA.markers,file = "Files/markers.Rdata")
save(scRNA, file = "Files/scRNA_acc_markers.Rdata")


##############Plot the proportion of different cancer cells in each cluster
load("Files/scRNA_acc_markers.Rdata")
table(scRNA$orig.ident)#Check the number of cells in each group
prop.table(table(Idents(scRNA)))
a<-table(scRNA$orig.ident,Idents(scRNA))#Number of cells in different groups
write.csv(a,file = "E:\\yjs\\bioinformation\\cancer\\Plots\\table\\cluster-sample.csv")
Cellratio <- prop.table(a, margin = 2)#Calculate the proportion of different cell populations in each group
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
table(Cellratio$Var2)
##Integrate GSM numbers by cancer cell type
BRAC<-Cellratio[c(1:5,17:21,33:37,49:53,65:69,81:85,97:101,113:117,129:133,145:149,161:165),]
LIHC<-Cellratio[c(6:7,22:23,38:39,54:55,70:71,86:87,102:103,118:119,134:135,150:151,166:167),]
LUAD<-Cellratio[c(8:16,24:32,40:48,56:64,72:80,88:96,104:112,120:128,136:144,152:160,168:176),]
table(BRAC$Var1)
table(LIHC$Var1)
table(LUAD$Var1)
FREQ_BRAC<- data.frame("0" = sum(BRAC[which(BRAC$Var2=="0"),3]),
                       "1" = sum(BRAC[which(BRAC$Var2=="1"),3]),
                       "2" = sum(BRAC[which(BRAC$Var2=="2"),3]),
                       "3" = sum(BRAC[which(BRAC$Var2=="3"),3]),
                       "4" = sum(BRAC[which(BRAC$Var2=="4"),3]),
                       "5" = sum(BRAC[which(BRAC$Var2=="5"),3]),
                       "6" = sum(BRAC[which(BRAC$Var2=="6"),3]),
                       "7" = sum(BRAC[which(BRAC$Var2=="7"),3]),
                       "8" = sum(BRAC[which(BRAC$Var2=="8"),3]),
                       "9" = sum(BRAC[which(BRAC$Var2=="9"),3]),
                       "10" = sum(BRAC[which(BRAC$Var2=="10"),3]))
FREQ_BRAC<-as.data.frame(t(FREQ_BRAC))
FREQ_BRAC <- data.frame(Var1 = "BRCA",
                        var2 = paste0("cluster", 0:10),
                        freq = FREQ_BRAC$V1)
FREQ_LIHC<- data.frame("0" = sum(LIHC[which(LIHC$Var2=="0"),3]),
                       "1" = sum(LIHC[which(LIHC$Var2=="1"),3]),
                       "2" = sum(LIHC[which(LIHC$Var2=="2"),3]),
                       "3" = sum(LIHC[which(LIHC$Var2=="3"),3]),
                       "4" = sum(LIHC[which(LIHC$Var2=="4"),3]),
                       "5" = sum(LIHC[which(LIHC$Var2=="5"),3]),
                       "6" = sum(LIHC[which(LIHC$Var2=="6"),3]),
                       "7" = sum(LIHC[which(LIHC$Var2=="7"),3]),
                       "8" = sum(LIHC[which(LIHC$Var2=="8"),3]),
                       "9" = sum(LIHC[which(LIHC$Var2=="9"),3]),
                       "10" = sum(LIHC[which(LIHC$Var2=="10"),3]))
FREQ_LIHC<-as.data.frame(t(FREQ_LIHC))
FREQ_LIHC <- data.frame(Var1 = "LIHC",
                        var2 = paste0("cluster", 0:10),
                        freq = FREQ_LIHC$V1)
FREQ_LUAD<- data.frame("0" = sum(LUAD[which(LUAD$Var2=="0"),3]),
                       "1" = sum(LUAD[which(LUAD$Var2=="1"),3]),
                       "2" = sum(LUAD[which(LUAD$Var2=="2"),3]),
                       "3" = sum(LUAD[which(LUAD$Var2=="3"),3]),
                       "4" = sum(LUAD[which(LUAD$Var2=="4"),3]),
                       "5" = sum(LUAD[which(LUAD$Var2=="5"),3]),
                       "6" = sum(LUAD[which(LUAD$Var2=="6"),3]),
                       "7" = sum(LUAD[which(LUAD$Var2=="7"),3]),
                       "8" = sum(LUAD[which(LUAD$Var2=="8"),3]),
                       "9" = sum(LUAD[which(LUAD$Var2=="9"),3]),
                       "10" = sum(LUAD[which(LUAD$Var2=="10"),3]))
FREQ_LUAD<-as.data.frame(t(FREQ_LUAD))
FREQ_LUAD <- data.frame(Var1 = "LUAD",
                        var2 = paste0("cluster", 0:10),
                        freq = FREQ_LUAD$V1)
freq<-rbind(FREQ_BRAC,FREQ_LIHC,FREQ_LUAD)
freq$var2 <- factor(freq$var2, levels = paste0("cluster", 0:10))#Set the factor level to prevent the ordinate from changing when drawing

library(ggplot2)
p<-ggplot(freq) + 
  geom_bar(aes(x =var2, y= freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cluster',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave("OUT/ratio_cluster.pdf",p,width = 10,height = 6)
ggsave("E:\\yjs\\bioinformation\\cancer\\Plots\\p1\\ratio_cluster.png",p,width = 10,height = 6)




##########################
###Mapping gene enrichment by marker
load("Files/markers.Rdata")
scRNA.markers<-scRNA.markers[which(scRNA.markers$p_val<0.05),]
library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(scRNA.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## Convert SYMBOL to ENTREZID
scRNA.markers=merge(scRNA.markers,ids,by.x='gene',by.y='SYMBOL')
scRNA.markers$cluster<-paste0("cluster",scRNA.markers$cluster)
View(scRNA.markers)
scRNA.markers$cluster<- factor(scRNA.markers$cluster, levels = paste0("cluster", 0:10))#Set the factor level to prevent the ordinate from changing when drawing

gcSample=split(scRNA.markers$ENTREZID, scRNA.markers$cluster) 

## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
#GO
aa<- compareCluster(gcSample,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db", ont= "BP",
                     pvalueCutoff=0.01,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05)
p <- dotplot(aa)
p1<- dotplot(aa,showCategory=5, includeAll=FALSE)
p1<- p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("OUT/KEGG_dotplot.pdf",p1,width = 9,height = 14)
ggsave("OUT/GO_dotplot.pdf",p1,width = 9,height = 32)
ggsave("E:\\yjs\\bioinformation\\cancer\\Plots\\p1\\GO_dotplot.png",p1,width = 8,height = 12)
ggsave("E:\\yjs\\bioinformation\\cancer\\Plots\\p1\\GO_dotplot.pdf",p1,width = 10,height =28)


#02resident tissue macrophages 13 TAM 68monocyte
load("Files/scRNA_acc_markers.Rdata")
library(Seurat)
library(dplyr)
library(tibble)
resident <- FindMarkers(object = scRNA, ident.1 = 0, ident.2 = 2, min.pct = 0.25)
resident <- resident %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))%>%rownames_to_column('SYMBOL')
TAM <- FindMarkers(object = scRNA, ident.1 = 1, ident.2 = 3, min.pct = 0.25)
TAM <- TAM %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% rownames_to_column('SYMBOL')
monocyto <- FindMarkers(object = scRNA, ident.1 = 6, ident.2 = 8, min.pct = 0.25)
monocyto <- monocyto%>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% rownames_to_column('SYMBOL')

library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

ENTREZID <- bitr(resident$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
resident <- resident %>% inner_join(ENTREZID,by= 'SYMBOL')

ENTREZID <- bitr(TAM$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
TAM <- TAM %>% inner_join(ENTREZID,by= 'SYMBOL')

ENTREZID <- bitr(monocyto$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
monocyto <- monocyto %>% inner_join(ENTREZID,by= 'SYMBOL')

genelist <- resident$avg_log2FC
names(genelist) <- resident$ENTREZID

kegg <- gseKEGG(genelist,
                organism = 'hsa',
                keyType = 'kegg',
                minGSSize = 10,
                pvalueCutoff = 0.05,
                seed = 123,pAdjustMethod = 'BH')
clus02<-data.frame(kegg)

library(enrichplot)
p <- ridgeplot(kegg,showCategory = 10,label_format = 60) 
#dir.create("diffgene&pathway")
write.csv(clus02,file = "diffgene&pathway/clus02.csv")
ggsave("diffgene&pathway/clus02.png",p,width = 8,height = 5)

genelist <- TAM$avg_log2FC
names(genelist) <- TAM$ENTREZID

kegg <- gseKEGG(genelist,
                organism = 'hsa',
                keyType = 'kegg',
                minGSSize = 10,
                pvalueCutoff = 0.05,
                seed = 123,pAdjustMethod = 'BH')
clus13<-data.frame(kegg)
p <- ridgeplot(kegg,showCategory = 10,label_format = 60) 
write.csv(clus13,file = "diffgene&pathway/clus13.csv")
ggsave("diffgene&pathway/clus13.png",p,width = 8,height = 5)
rm(kegg)

genelist <- monocyto$avg_log2FC
names(genelist) <- monocyto$ENTREZID

kegg <- gseKEGG(genelist,
                organism = 'hsa',
                keyType = 'kegg',
                minGSSize = 10,
                pvalueCutoff = 1,
                seed = 123,pAdjustMethod = 'BH')
clus68<-data.frame(kegg)
p <- ridgeplot(kegg,showCategory = 8,label_format = 60) 
#dir.create("diffgene&pathway")
write.csv(clus68,file = "diffgene&pathway/clus68.csv")
ggsave("diffgene&pathway/clus68.png",p,width = 8,height = 5)

###########
#GO
genelist <- resident$avg_log2FC
names(genelist) <- resident$SYMBOL

GO <- gseGO(genelist,ont = "BP",
              OrgDb = 'org.Hs.eg.db',
              keyType = 'SYMBOL',
                minGSSize = 10,
                pvalueCutoff = 0.05,
                seed = 123,pAdjustMethod = 'BH')

clus02<-data.frame(GO)
rm(GO)

genelist <- TAM$avg_log2FC
names(genelist) <- TAM$SYMBOL

GO <- gseGO(genelist,ont = "BP",
            OrgDb = 'org.Hs.eg.db',
            keyType = 'SYMBOL',
            minGSSize = 10,
            pvalueCutoff = 0.05,
            seed = 123,pAdjustMethod = 'BH')

clus13<-data.frame(GO)
rm(GO)

genelist <- monocyto$avg_log2FC
names(genelist) <- monocyto$SYMBOL

GO <- gseGO(genelist,ont = "BP",
            OrgDb = 'org.Hs.eg.db',
            keyType = 'SYMBOL',
            minGSSize = 10,
            pvalueCutoff = 1,
            seed = 123,pAdjustMethod = 'BH')

clus68<-data.frame(GO)

###########
#GSEA

library(stats)
genelist <- resident$avg_log2FC
names(genelist) <- resident$SYMBOL
hallmark <- clusterProfiler::read.gmt("macrophage/c7.all.v2023.2.Hs.entrez.gmt")
res <- GSEA(
  genelist,
  TERM2GENE = hallmark
)

clus02<-data.frame(res)
rm(res)

genelist <- TAM$avg_log2FC
names(genelist) <- TAM$SYMBOL

GO <- gseGO(genelist,ont = "BP",
            OrgDb = 'org.Hs.eg.db',
            keyType = 'SYMBOL',
            minGSSize = 10,
            pvalueCutoff = 0.05,
            seed = 123,pAdjustMethod = 'BH')

clus13<-data.frame(GO)
rm(GO)

genelist <- monocyto$avg_log2FC
names(genelist) <- monocyto$SYMBOL

res <- GSEA(
  genelist,
  TERM2GENE = hallmark
)

clus68<-data.frame(res)

###########
##cell cycle status
load("Files/scRNA_acc_markers.Rdata")
library(scran)
library(org.Hs.eg.db)
library(ggplot2)
library(readxl)
library(Seurat)
str(cc.genes)#seurat Package contains files

# Marker genes for the cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scRNA <- CellCycleScoring(scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#View cell cycle scores and phase assignments
head(scRNA[[]])

# Visualize the distribution of markers throughout the cell cycle
RidgePlot(scRNA, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

p<-DimPlot(scRNA,label = F,label.size = 3.0)
data<-scRNA@meta.data
data<-data[,c(7,16)]
cellcycle<-as.data.frame(table(data))
cellcycle<-data.frame(G1=cellcycle[1:11,3],
                 G2M=cellcycle[12:22,3],
                 S=cellcycle[23:33,3])
write.csv(cellcycle,file = "cellcycle.csv")
ggsave("OUT/UMAP_cellcycle.pdf",p,width = 5,height = 4)


##################
##Marker genes that produce proteins on the surface of cell membranes
####
library(readxl)
library(ggplot2)
data1<-read_excel("S2_File.xlsx",sheet = "Table A",col_names = TRUE)
load("Files/scRNA_acc_markers.Rdata")
load("Files/markers.Rdata")
genes<-data1$`ENTREZ gene symbol`
library(tidyverse)
top10 = scRNA.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
genes2<-top10$gene
gene<-genes[na.omit(match(genes2,genes))]

exp <- GetAssayData(scRNA, slot = "data")
head(scRNA$seurat_clusters)
cluster_info <- sort(scRNA$seurat_clusters)
head(cluster_info)

#####After splitting by group, the gene expression was averaged and then summarized
mat <- as.data.frame(exp[gene, names(cluster_info)])
mat0<-mat[,which(cluster_info=="0")]
mat1<-mat[,which(cluster_info=="1")]
mat2<-mat[,which(cluster_info=="2")]
mat3<-mat[,which(cluster_info=="3")]
mat4<-mat[,which(cluster_info=="4")]
mat5<-mat[,which(cluster_info=="5")]
mat6<-mat[,which(cluster_info=="6")]
mat7<-mat[,which(cluster_info=="7")]
mat8<-mat[,which(cluster_info=="8")]
mat9<-mat[,which(cluster_info=="9")]
mat10<-mat[,which(cluster_info=="10")]
##Confirm the total
sum(ncol(mat0),ncol(mat1),ncol(mat2),ncol(mat3),ncol(mat4),ncol(mat5),ncol(mat6),ncol(mat7),ncol(mat8),ncol(mat9),ncol(mat10))
sum(ncol(mat))
#########Take the average
mat0_mean = apply(mat0,1,mean)#1 Represents the average value of the display rows
mat1_mean = apply(mat1,1,mean)
mat2_mean = apply(mat2,1,mean)
mat3_mean = apply(mat3,1,mean)
mat4_mean = apply(mat4,1,mean)
mat5_mean = apply(mat5,1,mean)
mat6_mean = apply(mat6,1,mean)
mat7_mean = apply(mat7,1,mean)
mat8_mean = apply(mat8,1,mean)
mat9_mean = apply(mat9,1,mean)
mat10_mean = apply(mat10,1,mean)

mat_mean<- data.frame(Cluster0 = mat0_mean,
                      Cluster1 = mat1_mean,
                      Cluster2 = mat2_mean,
                      Cluster3 = mat3_mean,
                      Cluster4 = mat4_mean,
                      Cluster5 = mat5_mean,
                      Cluster6 = mat6_mean,
                      Cluster7 = mat7_mean,
                      Cluster8 = mat8_mean,
                      Cluster9 = mat9_mean,
                      Cluster10 = mat10_mean)


set.seed(123)
#Calculate the mean for each row
rowmean=apply(mat_mean,1,mean)
#Calculate the standard deviation for each row
rowsd=apply(mat_mean,1,sd)
#The gene expression value for each row is subtracted from the mean value of this row
data1=sweep(mat_mean,1,rowmean)
#The gene expression value for each row is divided by the standard deviation of this row
data2=sweep(data1,1,rowsd,'/')

bk = unique(c(seq(-1,1, length=100)))
library(pheatmap)
p<-pheatmap(data2,
         breaks = bk,
         show_colnames = T,
         angle_col = 90,
         border_color = NA,
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(colors = c("#0068D2","white","tomato"))(100),
         cellwidth = 10,cellheight = 10,
         fontsize=8)

ggsave("OUT/protein_heatmap.pdf",p,width = 5,height = 8)

##################
##Gene expression of the selected function
load("Files/scRNA_acc_markers.Rdata")
pro_inflammatory <- data.frame(gene = c("IL1B","CCL2","CCL3","CCL4","TNF","CXCL2","CXCL9","CXCL10","SOCS3"))
anti_inflammatory <- data.frame(gene = c("CD163","MRC1","MSR1","CCL13","CCL18","IL10","FOLR2","STAB1","CD163L1","SEPP1","F13A1","MERTK","AXL","CCL22","IL1RN","GPNMB","CHI3L1"))
MMP <- data.frame(gene = c("MMP1","MMP7","MMP9","MMP12","MMP14"))
others <- data.frame(gene = c("PPARG","FABP4","INHBA","ALDH2","LYVE1","EGFL7","CD209","CH25H","LILRB5"))

all<- rbind(pro_inflammatory,anti_inflammatory,MMP,others)


# Extract the original expression matrix
exp <- GetAssayData(scRNA, slot = "data")
head(scRNA$seurat_clusters)
cluster_info <- sort(scRNA$seurat_clusters)
head(cluster_info)

mat <- as.matrix(exp[, names(cluster_info)])
a<-match(rownames(mat),all$gene)
gene<-all[na.omit(a),]
mat <- as.data.frame(mat[gene, ])
b<-match(all$gene,rownames(mat))
mat<-mat[b,]
mat<- na.omit(mat)


mat0<-mat[,which(cluster_info=="0")]
mat1<-mat[,which(cluster_info=="1")]
mat2<-mat[,which(cluster_info=="2")]
mat3<-mat[,which(cluster_info=="3")]
mat4<-mat[,which(cluster_info=="4")]
mat5<-mat[,which(cluster_info=="5")]
mat6<-mat[,which(cluster_info=="6")]
mat7<-mat[,which(cluster_info=="7")]
mat8<-mat[,which(cluster_info=="8")]
mat9<-mat[,which(cluster_info=="9")]
mat10<-mat[,which(cluster_info=="10")]

sum(ncol(mat0),ncol(mat1),ncol(mat2),ncol(mat3),ncol(mat4),ncol(mat5),ncol(mat6),ncol(mat7),ncol(mat8),ncol(mat9),ncol(mat10))
sum(ncol(mat))

mat0_mean = apply(mat0,1,mean)#1Represents the average value of the display rows
mat1_mean = apply(mat1,1,mean)
mat2_mean = apply(mat2,1,mean)
mat3_mean = apply(mat3,1,mean)
mat4_mean = apply(mat4,1,mean)
mat5_mean = apply(mat5,1,mean)
mat6_mean = apply(mat6,1,mean)
mat7_mean = apply(mat7,1,mean)
mat8_mean = apply(mat8,1,mean)
mat9_mean = apply(mat9,1,mean)
mat10_mean = apply(mat10,1,mean)
mat_mean<- data.frame(Cluster0 = mat0_mean,
                      Cluster1 = mat1_mean,
                      Cluster2 = mat2_mean,
                      Cluster3 = mat3_mean,
                      Cluster4 = mat4_mean,
                      Cluster5 = mat5_mean,
                      Cluster6 = mat6_mean,
                      Cluster7 = mat7_mean,
                      Cluster8 = mat8_mean,
                      Cluster9 = mat9_mean,
                      Cluster10 = mat10_mean)


set.seed(123)
#Calculate the mean for each row
rowmean=apply(mat_mean,1,mean)
#Calculate the standard deviation for each row
rowsd=apply(mat_mean,1,sd)
#The gene expression value for each row is subtracted from the mean value of this row
data1=sweep(mat_mean,1,rowmean)
#The gene expression value for each row is divided by the standard deviation of this row
data2=sweep(data1,1,rowsd,'/')

bk = unique(c(seq(-1,1, length=100)))
library(pheatmap)
p<-pheatmap(data2,
         breaks = bk,
         show_colnames = T,
         angle_col = 90,
         border_color = NA,
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(colors = c("#0068D2","white","tomato"))(100),
         cellwidth = 10,cellheight = 10,
         fontsize=8,gaps_row = c(9,23,27))
ggsave("OUT/macro_gene_pheatmap.pdf",p,width = 4,height =7)

##############

#Alveolar macrophages kuppfer macrophages (HCC-specific macrophages) Macrophages expressed in different clusters within breast cancer
ARM <- c('FABP4','MCEMP1','MARCO')
KC <- c('Siglec1','Clec4f','Vsig4','Timd4')
FOLR <- c('SEPP1','MMP9','SLC40A1','CTSC', 'ADAMDEC1','TMEM176B','FOLR2','ENPP2','TSPAN4','PLD3','CXCL12','TMEM176A', 'MS4A4A','ATOX1',
          'FUCA1','MAF', 'IGF1','ITM2B','PLA2G2D','VCAM1','CREG1','TIMD4','STAB1','PDK4','PTGDS','CETP','IL18',
          'BLVRB','GPR34','CD163','F13A1','CCL18','LYVE1', 'LMNA','CXCL9','MS4A7','MRC1','CCL4'
)


install.packages('homologene')
library(homologene)

#Use the homologene for conversion
#@inTaxis the species number to which the entered gene list belongs, and 10090 is the mouse
#@outTaxis the species number to be converted into, and 9606 is human
KCh <- homologene(KC, inTax = 10090, outTax = 9606)[,2]
load("Files/scRNA_acc_markers.Rdata")
library(Seurat)
library(ggplot2)
a<-c(ARM,FOLR,KCh)
DefaultAssay(scRNA) <- "integrated"
p<-DotPlot(scRNA, features = unique(c(ARM,FOLR,KCh)))+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1 ,size = 8))


p <- DoHeatmap(scRNA, features = unique(c(ARM,FOLR,KCh)))+scale_fill_gradientn(colors = c("blue", "white", "red"))
ggsave("pheatmap_cancers_marker.png",p,width = 8,height = 7)




###############
#Differential expression of tissue-resident macrophages in Cluster0 Cluster2 among different tumors
setwd("E:/yjs/bioinformation/cancer/macrophage")
load("Files/scRNA_acc_markers.Rdata")
meta<-scRNA@meta.data
library(tibble)
library(dplyr)
meta<-rownames_to_column(meta)
cancer <- data.frame(orig.ident=c('GSM4476486','GSM4476487','GSM4476488','GSM4476489','GSM4476490','GSM5076749','GSM5076750','GSM5219674','GSM5219675','GSM5219676','GSM5219677','GSM5219678','GSM5219679','GSM5219680','GSM5219681','GSM5219682'),
                     cancertype = c(rep('BRCA',5),rep('LIHC',2),rep('LUAD',9))
                     )
meta<-right_join(meta,cancer,by = "orig.ident")
meta <- column_to_rownames(meta)
scRNA <- AddMetaData(scRNA,meta)
c0 = scRNA[,scRNA@meta.data$integrated_snn_res.0.5 %in% c(0)]
c2 = scRNA[,scRNA@meta.data$integrated_snn_res.0.5 %in% c(2)]
#cluster0 Differential genes for three different cancers
brca <- FindMarkers(c0, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "cancertype",
                     ident.1 ="BRCA")
lihc <- FindMarkers(c0, min.pct = 0.25, 
                    logfc.threshold = 0.25,
                    group.by = "cancertype",
                    ident.1 ="LIHC")
luad <- FindMarkers(c0, min.pct = 0.25, 
                    logfc.threshold = 0.25,
                    group.by = "cancertype",
                    ident.1 ="LUAD")
brca <- brca %>% filter(p_val < 0.05) %>% arrange(desc(avg_log2FC))
lihc <- lihc %>% filter(p_val < 0.05) %>% arrange(desc(avg_log2FC))
luad <- luad %>% filter(p_val < 0.05) %>% arrange(desc(avg_log2FC))

#cluster2 Differential genes for three different cancers
brca2 <- FindMarkers(c2, min.pct = 0.25, 
                    logfc.threshold = 0.25,
                    group.by = "cancertype",
                    ident.1 ="BRCA")
lihc2 <- FindMarkers(c2, min.pct = 0.25, 
                    logfc.threshold = 0.25,
                    group.by = "cancertype",
                    ident.1 ="LIHC")
luad2 <- FindMarkers(c2, min.pct = 0.25, 
                    logfc.threshold = 0.25,
                    group.by = "cancertype",
                    ident.1 ="LUAD")
brca2 <- brca2 %>% filter(p_val < 0.05) %>% arrange(desc(avg_log2FC))
lihc2 <- lihc2 %>% filter(p_val < 0.05) %>% arrange(desc(avg_log2FC))
luad2 <- luad2 %>% filter(p_val < 0.05) %>% arrange(desc(avg_log2FC))


############
library('devtools')
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
brca<-rownames_to_column(brca)
brca$cluster <- "BRCA"
lihc<-rownames_to_column(lihc)
lihc$cluster <- "LIHC"
luad<-rownames_to_column(luad)
luad$cluster <- "LUAD"
diff<-rbind(brca,lihc,luad)
colnames(diff)[1]<-"gene"
head(diff)
jjVolcano(diffData = diff,
          log2FC.cutoff = 0.25, 
          col.type = "adjustP",#The graph color is mapped to the up-down type, while adjustP is mapped to whether the corrected P value is less than 0.01
          size  = 3.5, #Sets the size of the point
          #myMarkers = c("CD68","CD14"),#Provide a gene name and set your own gene
          fontface = 'italic', #Set the font style
          #aesCol = c('purple','orange'), #The color of the set
          tile.col = corrplot::COL2('RdBu', 15)[c(6,9,11)], #Set the color of the cluster
          #col.type = "adjustP", #Set the correction method
          topGeneN = 5,#Set up the gene that exhibits topN
          polar = F, #Adjust the shape
          legend.position = c(0.82,0.95)
)

brca2<-rownames_to_column(brca2)
brca2$cluster <- "BRCA"
lihc2<-rownames_to_column(lihc2)
lihc2$cluster <- "LIHC"
luad2<-rownames_to_column(luad2)
luad2$cluster <- "LUAD"
diff<-rbind(brca2,lihc2,luad2)
colnames(diff)[1]<-"gene"
head(diff)
jjVolcano(diffData = diff,
          log2FC.cutoff = 0.25, 
          col.type = "adjustP",#The graph color is mapped to the up-down type, while adjustP is mapped to whether the corrected P value is less than 0.01
          size  = 3.5, #Sets the size of the point
          #myMarkers = c("CD68","CD14"),#Provide a gene name and set your own gene
          fontface = 'italic', #Set the font style
          #aesCol = c('purple','orange'), #Sets the color of the point
          tile.col = corrplot::COL2('RdBu', 15)[c(6,9,11)], #Set the color of the cluster
          #col.type = "adjustP", #Set the correction method
          topGeneN = 5,#Set up the gene that exhibits topN
          polar = F, #Adjust the shape
          legend.position = c(0.82,0.98)
)

