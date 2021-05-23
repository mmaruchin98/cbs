library(Rtsne)
library(data.table)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)


pbmc.data <- Read10X(data.dir = "/Users/Micha³/Desktop/filtered_matrices_mex/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc68k", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc.embedet <- Embeddings(object = pbmc, reduction = "pca")[,1:20]
k.mean <- kmeans(pbmc.embedet, 10)
k.mean.clusters <- k.mean$cluster
set.seed(1)
tsne_out <- Rtsne(pbmc.embedet,pca=F,perplexity=30)
tsne_out_pos = data.table(tsne_out$Y)
tsne_out_pos$cluster <- k.mean$cluster
ggplot(tsne_out_pos) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))  + labs(color = "cluster")

c1 <- pbmc.embedet[k.mean.clusters == 1,]
c2 <- pbmc.embedet[k.mean.clusters == 2,]
c3 <- pbmc.embedet[k.mean.clusters == 3,]
c4 <- pbmc.embedet[k.mean.clusters == 4,]
c5 <- pbmc.embedet[k.mean.clusters == 5,]
c6 <- pbmc.embedet[k.mean.clusters == 6,]
c7 <- pbmc.embedet[k.mean.clusters == 7,]
c8 <- pbmc.embedet[k.mean.clusters == 8,]
c9 <- pbmc.embedet[k.mean.clusters == 9,]
c10 <- pbmc.embedet[k.mean.clusters == 10,]

png(file="lab8_c1.png")
c1.kmean <- kmeans(c1, 5)
clusters.c1 <- c1.kmean$cluster
c1.tsne <- Rtsne(c1, pca = F, perplexity = 30)
c1.tsne.dt <- data.table(c1.tsne$Y) 
c1.tsne.dt$cluster <- c1.kmean$cluster
ggplot(c1.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()

png(file="lab8_c2.png")
c2.kmean <- kmeans(c2, 5)
clusters.c2 <- c2.kmean$cluster
c2.tsne <- Rtsne(c2, pca = F, perplexity = 30)
c2.tsne.dt <- data.table(c2.tsne$Y) 
c2.tsne.dt$cluster <- c2.kmean$cluster
ggplot(c2.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()
 
png(file="lab8_c3.png")
c3.kmean <- kmeans(c3, 5)
clusters.c3 <- c3.kmean$cluster
c3.tsne <- Rtsne(c3, pca = F, perplexity = 30)
c3.tsne.dt <- data.table(c3.tsne$Y) 
c3.tsne.dt$cluster <- c3.kmean$cluster
ggplot(c3.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()

png(file="lab8_c4.png")
c4.kmean <- kmeans(c4, 5)
clusters.c4 <- c4.kmean$cluster
c4.tsne <- Rtsne(c4, pca = F, perplexity = 30)
c4.tsne.dt <- data.table(c4.tsne$Y) 
c4.tsne.dt$cluster <- c4.kmean$cluster
ggplot(c4.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()


png(file="lab8_c5.png")
c5.kmean <- kmeans(c5, 5)
clusters.c5 <- c5.kmean$cluster
c5.tsne <- Rtsne(c5, pca = F, perplexity = 30)
c5.tsne.dt <- data.table(c5.tsne$Y) 
c5.tsne.dt$cluster <- c5.kmean$cluster
ggplot(c5.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()

png(file="lab8_c6.png")
c6.kmean <- kmeans(c6, 5)
clusters.c6 <- c6.kmean$cluster
c6.tsne <- Rtsne(c6, pca = F, perplexity = 30)
c6.tsne.dt <- data.table(c6.tsne$Y) 
c6.tsne.dt$cluster <- c6.kmean$cluster
ggplot(c6.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()

png(file="lab8_c7.png")
c7.kmean <- kmeans(c7, 5)
clusters.c7 <- c7.kmean$cluster
c7.tsne <- Rtsne(c7, pca = F, perplexity = 30)
c7.tsne.dt <- data.table(c7.tsne$Y) 
c7.tsne.dt$cluster <- c7.kmean$cluster
ggplot(c7.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()

png(file="lab8_c8.png")
c8.kmean <- kmeans(cluster.8, 5)
clusters.cluster.8 <- cluster.8.kmean$cluster
cluster.8.tsne <- Rtsne(cluster.8, pca = F, perplexity = 30)
cluster.8.tsne.dt <- data.table(cluster.8.tsne$Y) 
cluster.8.tsne.dt$cluster <- cluster.8.kmean$cluster
ggplot(cluster.8.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()

png(file="lab8_c9.png")
c9.kmean <- kmeans(c9, 5)
clusters.c9 <- c9.kmean$cluster
c9.tsne <- Rtsne(c9, pca = F, perplexity = 30)
c9.tsne.dt <- data.table(c9.tsne$Y) 
c9.tsne.dt$cluster <- c9.kmean$cluster
ggplot(c9.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()

png(file="lab8_c10.png")
c10.kmean <- kmeans(c10, 5)
clusters.c10 <- c10.kmean$cluster
c10.tsne <- Rtsne(c10, pca = F, perplexity = 30)
c10.tsne.dt <- data.table(c10.tsne$Y) 
c10.tsne.dt$cluster <- c10.kmean$cluster
ggplot(c10.tsne.dt) + geom_point(aes(x=V1, y=V2, col = as.factor(cluster)))
dev.off()


