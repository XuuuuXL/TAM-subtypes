library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)

load("macrophage/Files/scRNA_acc_markers.Rdata")
data.input  <- scRNA@assays$RNA@data
identity = data.frame(group =factor(paste0("clus",scRNA$seurat_clusters),levels = paste0("clus",0:10)) , row.names = names(scRNA$seurat_clusters)) # create a dataframe consisting of the cell labels

unique(identity$group) # check the cell labels
#Create a Cell Chat object and enter the data after log
cellchat <- createCellChat(data.input)
cellchat

str(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
#[1] 2738 1300 1299  466  446  425  408  252  208  182   77


CellChatDB <- CellChatDB.human#CellChatDB.mouse

str(CellChatDB) #View database information
#It contains four dataframes: interaction, complex, cofactor, and geneInfo
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
#After the ligand receptor relationship is reached, projectData projects the expression value of the ligand receptor pair onto the PPI to correct the expression value in @data.signaling. The results are saved in @data.project


cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")


#Aggregate communication networks between cells can be calculated by counting the number of links or by summarizing communication probabilities
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")

######

#Calculate and summarize the sum of intensities in different pathways
df <- read.csv('net_pathway.csv',header = T,row.names = 1)
dfsum <- aggregate(df$prob, by=list(type=df$pathway_name),sum)

#Draw a pie chart
#colour
library(randomcoloR)
set.seed(123)
palette <- distinctColorPalette(58) #60 species with significant differences
#Build data
dfsum$ratio <- dfsum$x/sum(df$prob)
dfsum <- dfsum[order(dfsum$x,decreasing = T),]
a <- data.frame(type = dfsum$type,
                ratio = dfsum$ratio)


pie(dfsum$ratio, labels=dfsum$type,
    radius = 1.0,clockwise=T,
    main = " ")

library(ggplot2)
library(ggforce)
index=order(a$ratio)


ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(),
        legend.key.size = unit(35,'pt'),
        legend.margin = margin(2,2,2,2,'cm'),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+scale_fill_manual(values = palette)+
  geom_arc_bar(data=a,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type)
  )


cellchat <- aggregateNet(cellchat)

cellchat@netP$pathways
head(cellchat@LR$LRsig)
cellchat@netP$pathways
levels(cellchat@idents) 
save(cellchat,file = "cellchat.rdata")
#Count how many of each type of cell there are
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
p1<-netVisual_circle((cellchat@net$count), vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
p2<-netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
ggsave("net_number.png",plot =  p1, width = 7, height = 7,device = "png")
ggsave("net_strength.png",plot =  p2, width = 7, height = 7)


#Check the signals emitted by each type of cell
mat <- cellchat@net$count
par(mfrow = c(3,4), xpd=TRUE)
p<-list()
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  p[[i]]<-netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
 
  }
p[[11]]

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=T)
p<-list()
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  p[[i]]<-netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
p[[11]]


#Visualization of single signaling pathways or ligand-receptor-mediated cellular interactions (hierarchical graphs, network diagrams, chord diagrams, heat maps)
cellchat@netP$pathways  #See what signaling pathways are available
# [1] "SPP1"       "GALECTIN"   "CCL"        "MIF"        "VISFATIN"   "IL1"       
# [7] "CXCL"       "BAFF"       "CSF"        "RESISTIN"   "VEGF"       "TGFb"      
# [13] "GAS"        "IFN-II"     "GRN"        "COMPLEMENT" "TNF"        "IL4"       
# [19] "LT"         "ANNEXIN"    "IL2"        "IL10"       "SEMA3"      "BTLA"       
# [25] "LIGHT"      "IL6"        "IGF"        "GDF"        "MK"         "UGRP1"     
# [31] "NTS"        "EGF"        "ACTIVIN"    "PROS"       "CSF3"       "CD40"      
# [37] "LIFR"       "TWEAK"      "OSM"        "BMP"        "TRAIL"      "PSAP"      
# [43] "CD137"      "FASLG"      "OX40"       "NRG"        "PRL"        "GH"        
# [49] "HGF"        "EPO"        "FGF"        "KIT"        "NT"         "PARs"      
# [55] "NGF"        "OPIOID"     "FLT3"       "WNT"

pathways.show <- c("SPP1")  

#Hierarchy chart
levels(cellchat@idents)    # show all celltype
# [1] [1] "clus0"  "clus1"  "clus2"  "clus3"  "clus4"  "clus5"  "clus6"  "clus7"  "clus8"  "clus9"  "clus10"   
vertex.receiver = c(1,2,3) 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize,layout = "hierarchy")


#Network diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# save as TIL/SPP1_circle.pdf
ggsave("SPP1_circle.png", p, width = 10, height = 10)


#Calculate the contribution of the ligand acceptor to the selected signaling pathway (in this case, see which signaling pathway contributes the most to SPP1)
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.SPP1 <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE) #Extract all ligand receptors that contribute to SPP1

#Extract the ligand-receptor pairs that contribute the most to this pathway for presentation (other ligand receptor pairs can also be selected)
LR.show <- pairLR.SPP1[1,] 

vertex.receiver = c(1,2,3) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout = "hierarchy")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show,layout = "chord")


#Automatically save the interaction results of each signaling path in batches
# Access all the signaling pathways showing significant communications. Find out all the signaling pathways
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(1,2,3) 
dir.create("all_pathways_com_circle") #Create a folder to save the batch drawing results
setwd("all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
            vertex.receiver = vertex.receiver, layout = "circle") #Draw a network diagram
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
         plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
}
setwd("../")

#Visualization of multiple ligand-receptor-mediated cellular interactions
#All ligand receptors
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#Recipient cells and ligand cells need to be specified
p = netVisual_bubble(cellchat, sources.use = c(1,2,3), 
                     targets.use = c(4,5,6,7,8,9,10,11), remove.isolate = FALSE)
ggsave("TIL/Cluster_bubble.pdf", p, width = 10, height = 18)
ggsave("E:\\yjs\\bioinformation\\cancer\\Plots\\p2\\Cluster_bubble.png", p, width = 10, height = 18)
# save as TIL/cluster_bubble.pdf

p = netVisual_bubble(cellchat, sources.use = c(4,5,6,7,8,9,10,11), 
                     targets.use = c(1,2,3), remove.isolate = FALSE)
ggsave("TIL/Cluster_bubble2.pdf", p, width = 10, height = 18) 

# save as TIL/cluster_bubble.pdf


#Specify the ligand receptor
p<-netVisual_bubble(cellchat, sources.use = c(1,2,3), targets.use = c(5,6,7,8,9,10,11), 
                 signaling = c("SPP1","GALECTIN","MIF"), remove.isolate = FALSE)
ggsave("Cluster_bubble_signal.png", p, width = 10, height = 10)

p<-netVisual_bubble(cellchat, sources.use = c(5,6,7,8,9,10,11), targets.use =c(1,2,3) , 
                 signaling = c("SPP1","GALECTIN","MIF"), remove.isolate = FALSE)
ggsave("Cluster_bubble_signal2.png", p, width = 10, height = 10)

#Specify signaling pathways or ligand-receptors and specify cells
# show all the significant interactions (L-R pairs) based on user's input
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("SPP1","GALECTIN","MIF"))
netVisual_bubble(cellchat, sources.use = c(1,2,3), targets.use = c(5,6,7,8,9,10,11), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)
#Expression of all genes involved in a signaling pathway in a cell population
## Plot the signaling gene expression distribution
p = plotGeneExpression(cellchat, signaling = "SPP1")
ggsave("SPP1_GeneExpression_vln.pdf", p, width = 8, height = 8)
ggsave("SPP1_GeneExpression_vln.png", p, width = 8, height = 8)

p = plotGeneExpression(cellchat, signaling = "SPP1", type = "dot",color.use = "rainbow")
ggsave("SPP1_GeneExpression_dot.pdf", p, width = 8, height = 6)

#Social Network Analysis (Role Recognition in Communication Networks)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#Identify the role/role of each cell class in the signaling pathway by calculating the network centrality metrics for each cell population
netAnalysis_signalingRole_network(cellchat, signaling = c("SPP1","GALECTIN","MIF"), 
                                  width = 15, height = 6, font.size = 10)


#Non-negative matrix factorization (NMF) to identify cellular communication patterns (here is a more standard application of NMF)
#Pattern recognition of signal output cells
#It is appropriate to break down the calculation into several patterns (this step is slow). When using NMF to subpopulation cells for subpopulation segmentation, it is better to choose a value that is a little more than the cell type if not tested)
library(NMF)
library(ggalluvial)
p<-selectK(cellchat, pattern = "outgoing")
ggsave("NMF_outgoing_selectK.png", p, width = 10, height = 5)

nPatterns = 3 # Pick the first point in the curve that goes down (it starts to fall from 2)
dev.off()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, 
                                          width = 4, height = 11, font.size = 6)

#Alluvial map
p<-netAnalysis_river(cellchat, pattern = "outgoing")
ggsave("NMF_outgoing_comPattern_river.png", p, width = 10, height = 9)
netAnalysis_dot(cellchat, pattern = "outgoing")
#Pattern recognition of signal input cells
p<-selectK(cellchat, pattern = "incoming") 
ggsave("NMF_incoming_selectK.png", p, width = 10, height =5)
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, 
                                          width = 5, height = 9, font.size = 6)
p<-netAnalysis_river(cellchat, pattern = "incoming")
ggsave("NMF_incoming_comPattern_river.png", p, width = 10, height = 9)

netAnalysis_dot(cellchat, pattern = "incoming")

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#Error in do_one(nmeth) : NA/NaN/Inf in foreign function call (arg 1)
p = netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
ggsave("Manifold_functional_cluster.pdf", p, width = 8, height = 6)


#Epidemiology and classification based on topological similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
p = netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
ggsave("Manifold_structural_cluster.pdf", p, width = 8, height = 6)

## The end
saveRDS(cellchat, file = "cellchat.rds")



###################################

library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(devtools)
#Use the merge function to merge seurat objects
scRNAlist <- list()

scRNAlist[[1]] <- scRNA_BRCA
scRNAlist[[2]] <- scRNA_LIHC
scRNAlist[[3]] <- scRNA_LUAD

names(scRNAlist)[1]<-"BRCA"
names(scRNAlist)[2]<-"LIHC"
names(scRNAlist)[3]<-"LUAD"
#scRNAlist is a list of seurat objects saved by the previous code run
#Data integration is preceded by data normalization and selection of high-variable genes for each sample's seurat object

for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}
##Finding anchors based on VariableFeatures takes a long time
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
##Consolidate data with anchors and take longer to run
scRNA <- IntegrateData(anchorset = scRNA.anchors)
DefaultAssay(scRNA) <- "integrated"
dim(scRNA)
#[1] 2000 57766

save(scRNA,file = "Files/scRNA_acc.Rdata")
library(ggplot2)
#Dimensionality reduction clustering
scRNA1 <- ScaleData(scRNA, features = VariableFeatures(scRNA))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1),seed.use = 1223)
plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA1, ndims=50, reduction="pca") 
plotc <- plot1+plot2
#dir.create("OUT")
ggsave("OUT/pca.pdf", plot = plotc, width = 8, height = 4)
#Select the principal component
pc.num=1:30


scRNA1 <- FindNeighbors(scRNA1, dims = pc.num)
scRNA_cluster <- scRNA1
library(clustree)
scRNA_cluster <- FindClusters(
  object = scRNA_cluster,
  resolution = c(seq(.2,1.2,.1))
)
p<-clustree(scRNA_cluster@meta.data, prefix = "integrated_snn_res.")

scRNA1 <- FindClusters(scRNA1, resolution = 0.7)
head(Idents(scRNA1),5)

table(scRNA1@meta.data$seurat_clusters)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 7257 5080 4567 4084 3907 2807 2288 2156 2149 1928 1805 1801 1554 1339 1283 1272 1269 1097  979 
# 19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34 
# 970  945  895  846  720  716  674  594  580  474  455  409  388  213  138  127 
metadata <- scRNA1@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)


#tSNE
scRNA1 = RunTSNE(scRNA1, dims = pc.num,seed.use = 1223)
embed_tsne <- Embeddings(scRNA1, 'tsne')   #Extract TSNE plot coordinates
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
scRNA1 <- RunUMAP(scRNA1, dims = pc.num,seed.use = 1223)
embed_umap <- Embeddings(scRNA1, 'umap')   #Extract the coordinates of the UMAP diagram
#write.csv(embed_umap,'cluster1/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
ggsave("OUT/UMAP_cluster.pdf", plot = plot3, width = 8, height = 7)

#group_by_sample
plot4 = DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')
ggsave("OUT/UMAP_sample.pdf", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("OUT/UMAP_cluster_sample.pdf", plot = plotc, width = 10, height = 5)


library(patchwork)
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("OUT/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
save(scRNA1,file = "Files/scRNA_acc_umap.Rdata")


####
library(SingleR)
refdata <- get(load("HumanPrimaryCellAtlasData.Rdata"))
assay(refdata)[1:4,1:4]
head(refdata@colData)
head(refdata)
library(celldex)
library(Seurat)
DefaultAssay(scRNA1) <- "RNA"
testdata <- GetAssayData(scRNA1, slot="data")
clusters <- scRNA1@meta.data$integrated_snn_res.0.7
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts"
)

table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) #The following is the identification result of singleR cell cluster.
scRNA1@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA1@meta.data[which(scRNA1@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
  }
p1 <- DimPlot(scRNA1, group.by="celltype", label=F )
p1
ggsave("OUT/celltype_UMAP.png",p1,width = 8,height = 7)
ggsave("OUT/celltype_UMAP.pdf",p1,width = 8,height = 7)
save(scRNA1,file = "Files/scRNA_celltype.Rdata")


##########CELLCHAT
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
load("Files/scRNA_celltype.Rdata") 
data.input = scRNA1@assays$RNA@data # normalized data matrix

identity = data.frame(group =scRNA1$celltype, row.names = names(scRNA1$celltype))
head(identity)
unique(identity$group) # check the cell labels
cellchat <- createCellChat(data.input)
cellchat
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
#[1]  3101   388   594  8898  1272  8331  5346   970  3854   409 23320  1283
CellChatDB <- CellChatDB.human

str(CellChatDB) #View database information
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
showDatabaseCategory(CellChatDB)

#In CellChat, you can select specific information to describe cell-cell interactions, and characterize cell-cell interactions from specific aspects
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)

#Pre-process expression data for cell communication analysis

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 5) # do parallel

#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)


#Inference of cellular communication networks

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "Files/net_lr.csv")

#Inference of cell-cell communication at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

df.netp <- subsetCommunication(cellchat, slot.name = "netP")
df <- df.netp[df.netp$source=='Macrophage'&df.netp$target=='Epithelial_cells',]
df2 <- df.netp[df.netp$source=='Epithelial_cells'&df.netp$target=='Macrophage',]
df <- rbind(df,df2)
#Computationally integrated cellular communication networks
cellchat <- aggregateNet(cellchat)
save(cellchat,file = 'Files/scRNA_cellchat.Rdata')

cellchat@netP$pathways
head(cellchat@LR$LRsig)
cellchat@netP$pathways
levels(cellchat@idents) 

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
#Number of interactions or total interaction intensity (specific gravity)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Visualization of cellular communication networks

pathways.show <- c("SPP1") 
vertex.receiver = seq(1,4) 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
p<-netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
ggsave( 'SPP1pathway_chord.png',plot = p,width = 8, height = 7)

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#
#Calculate the contribution of each ligand receptor to the overall signaling pathway and visualize cellular communication regulated by a single ligand receptor pair
p<-netAnalysis_contribution(cellchat, signaling = pathways.show)
ggsave( 'SPP1pathway_L-R.png',plot = p,width = 8, height = 7)
#Cell-to-cell communication regulated by a single ligand receptor pair
pairLR.SPP1 <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.SPP1[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
p <-netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
ggsave( 'SPP1-CD44.png',plot = p,width = 8, height = 7)

p<-plotGeneExpression(cellchat, signaling = "SPP1")
ggsave( 'SPP1vln.png',plot = p,width = 8, height = 7)




###################
library(TCGAbiolinks)
library(tidyr)
projects <- getGDCprojects()
library(dplyr)
projects <- projects %>% 
  as.data.frame() %>% 
  dplyr::select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))

CD44<-list()
ITGB1 <-list()
ITGB6 <- list()
for (i in 1:33) {
  filename1 = paste0("E:/yjs/bioinformation/cancer/KM/result/",projects$project_id[i],"-tpm.txt")
  mRNA<-read.table(filename1,header = T,check.names = F)
  mRNA<-mRNA[,-c(1,3)]
  data <- aggregate(x = mRNA[,-1],
                    by = list(mRNA$gene_name),
                    FUN = mean)
  # Log conversion of TPM data
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
  a <- as.data.frame(t(logTPM))
b<-which(rownames(a)== 'gene_name')
colnames(a) <- a['gene_name',]
a<-a[-b,]
y <- as.numeric(a[,'CD44'])
colnames <- colnames(a)
cor_data_df1 <- data.frame(colnames)
for (j in 1:length(colnames)){
  test <- cor.test(as.numeric(a[,i]),y,type="spearman")
  cor_data_df1[j,2] <- test$estimate
  cor_data_df1[j,3] <- test$p.value
}
names(cor_data_df1) <- c("symbol","correlation","pvalue")
CD44[[i]] <- cor_data_df1

y <- as.numeric(a[,'ITGB1'])
colnames <- colnames(a)
cor_data_df2 <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(a[,i]),y,type="spearman")
  cor_data_df2[i,2] <- test$estimate
  cor_data_df2[i,3] <- test$p.value
}
names(cor_data_df2) <- c("symbol","correlation","pvalue")
ITGB1[[i]] <- cor_data_df2

y <- as.numeric(a[,'ITGB6'])
colnames <- colnames(a)
cor_data_df3 <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(a[,i]),y,type="spearman")
  cor_data_df3[i,2] <- test$estimate
  cor_data_df3[i,3] <- test$p.value
}
names(cor_data_df3) <- c("symbol","correlation","pvalue")
ITGB6[[i]] <- cor_data_df3
}

