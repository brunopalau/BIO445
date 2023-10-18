##Starterscript  

##################################################
##### Installation of libraries ##################
##################################################

##Install BiocManager
##For the prompts use "a" or "yes" (whatever fits)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#Part 2
##Install Packages
##For the prompts use "a" or "yes" (whatever fits)
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("edgeR")
BiocManager::install("AnnotationDbi")
#BiocManager::install(organism, character.only = TRUE)
install.packages(ggridges)

#load libraries
#library(organism, character.only = TRUE)
library(clusterProfiler)
library("org.Hs.eg.db")
library(enrichplot)
library(GenomeInfoDb)
library(org.Hs.eg.db)
library(TBSignatureProfiler)
library(SummarizedExperiment)
library(pathview)
library(ggridges)
library(AnnotationDbi)
library(edgeR)

#Part 1
##Install Packages
##For the prompts use "a" or "yes" (whatever fits)
BiocManager::install("TBSignatureProfiler")
install.packages(tidyverse)
install.packages(matrixStats)
install.packages(limma)
#Load libraries
library(matrixStats)
library(tidyverse)
library(limma)

###################################################
##### Part 1: RNAseq data manual ##################
###################################################

##Load raw data data
hivtb_data = TB_hiv

##Extract gene names and counts
counts = hivtb_data@assays@data@listData[["counts"]]

##################################################
##### Your code or from the script ###############
##################################################



library(ggplot2)
library(tidyverse)
library(tidyr)


logcounts_2 <- log(counts)

boxplot(logcounts_2)

probemeans_2 <- rowMeans(logcounts_2)

q25_2 <- quantile(probemeans_2, 0.25)
whichtosave_2 <- which(probemeans_2 > q25)
q25logdata_2 <- logcounts[whichtosave_2,]

mydata_2 = q25logdata_2[apply(q25logdata_2,1,IQR) > 1,]

# do pca
tdata_2 <- t(mydata_2)

pca_2 = prcomp(tdata_2,scale=T)


# create dataframe
df_2 <- as.data.frame(pca_2$x)
conditions = c(rep("TB",16),rep("NOTB",15))

clusters_2 <- kmeans(df_2[,1:31],2)$cluster

df_2$condition = conditions
df_2$clusters = as.factor(clusters_2)


ggplot(df_2, aes(x = PC1,y = PC2, color = clusters, shape = condition))+
  geom_point()+
  ggtitle("pca on log2")
ggsave("log_2_pca.png")

ggplot(df_2, aes(x = condition, fill = clusters))+
  geom_bar()+
  ggtitle("clusters with log2")
ggsave("log_2_clusters.png")


# log 10
logcounts <- log10(counts)

boxplot(logcounts)

probemeans <- rowMeans(logcounts)

q25 <- quantile(probemeans, 0.25)
whichtosave <- which(probemeans > q25)
q25logdata <- logcounts[whichtosave,]

mydata = q25logdata[apply(q25logdata,1,IQR) > 1,]

# do pca
tdata <- t(mydata)

pca = prcomp(tdata,scale=T)

# create dataframe
df <- as.data.frame(pca$x)
conditions = c(rep("TB",16),rep("NOTB",15))

df$condition = conditions

clusters <- kmeans(df[,1:31],2)$cluster
df$clusters = as.factor(clusters)


ggplot(df, aes(x = PC1,y = PC2, color = clusters, shape = condition))+
  geom_point()+
  ggtitle("pca on log10")
ggsave("log_10_pca.png")

ggplot(df, aes(x = condition, fill = clusters))+
  geom_bar()+
  ggtitle("clusters with log10")
ggsave("log_10_clusters.png")



# do clustering on pca
library(tidyr)
library(tidyr)
library(tidyverse)
clusters_2 <- kmeans(df%>%
                       select(c(-condition)),2)









library(# pearson correlation
pearsonCorr = as.dist(1 - cor(mydata))
hC = hclust(pearsonCorr)
plot(hC, labels = conditions)

# euclidean distance

euclidean <- function(a, b){return (sqrt(sum((a - b)^2)))}
mat <- matrix(nrow=nrow(pca$x),ncol = nrow(pca$x))


for (i in 1:nrow(pca$x)){
  row <- pca$x[i,]
  for (j in 1:nrow(pca$x)){
    dist = euclidean(row,pca$x[j,])
    mat[i,j] <- dist
  }
}

euc_dist = as.dist(mat)
hC = hclust(euc_dist)
plot(hC, labels = conditions)


library(ggraph)
library(igraph)
library(tidyverse)

# Create a graph object 
mygraph <- graph_from_data_frame(df)

# Basic tree
ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal() +
  geom_node_point() +
  theme_void()

heatmap(mydata)

################################################################################
##### Part 2: Analysis with predefined pipelines from R packages ###############
################################################################################

##Alternative analysis with edgeR
##Load data again
hivtb_data = TB_hiv
##Extract counts per million 
hivtb_data = mkAssay(hivtb_data, log = TRUE, counts_to_CPM = TRUE)

##Look for genes with at least 3 samples with counts above 0.5 count per million
cpmdat = rownames(hivtb_data)[rowSums(assays(hivtb_data)$cpm > 0.5) >= 3]

##extract raw counts for each genes and filter those out with less than 3 genes with 0.5 counts per million, 
##transform all samples columns into numeric variables
countdat = assays(hivtb_data)$counts[rowSums(assays(hivtb_data)$cpm > 0.5) >= 3,]

##create edgeR object (group "1" = TB, group "0" = noTB)
egr = DGEList(counts = as.matrix(countdat[,1:31]), group = c(rep(1,16),rep(0,15)))

##normalize fold changes in gene expression between samples
egr = calcNormFactors(egr)

##Dispersion / distribution of gene counts
egr = estimateDisp(egr)

##Common dispersion is assuming equal mean and standard-deviation among genes
##Tagwise dispersion assumes gene-specific dispersions (or of genes with a similar distribution)
plotBCV(egr)

##Test gene expression differences between TB and noTB 
##The test is based on a negative binomial distribution 
et = exactTest(egr)

et$table %>% 
  mutate(logp = -log10(PValue)) %>%
  ggplot(., aes(logFC,logp )) +
  geom_point(aes(color = logp), size = 2/5) +
  xlab(expression("log2 fold change")) + 
  ylab(expression("-log10 pvalue")) +
  scale_color_viridis_c()


####################################################
########### Gene enrichment analysis ###############
####################################################

#We need to specify our organism of interest, so humans.
organism = "org.Hs.eg.db"

##we need log2foldchange from the previous analysis and gene names
genes = cbind(et$table$logFC,cpmdat$V1) %>% as.data.frame()

#The gene names, as they are right now are in the "Gene card symbol" format.
#For the analysis we should changed them to the ensembl coding.
hs = org.Hs.eg.db
genenames = select(hs, 
                   keys = genes$V2,
                   columns = c("ENSEMBL", "SYMBOL"),
                   keytype = "SYMBOL",
                   multiVals = "First")

#Some genes could not be recognized and are NA. Other have multiple ensembl ids.
#For simplicity remove the NAs and duplicates (take the first).
genes = genes %>% merge(.,genenames, by.x = "V2", by.y = "SYMBOL") %>% filter(ENSEMBL != "<NA>") %>% filter(!duplicated(ENSEMBL))

#Now extract the log2fold changes into a single vector, 
#name each value with the corresponding ensembl gene name, and sort the values decreasingly.
genenrich = as.numeric(genes$V1)
names(genenrich) = genes$ENSEMBL

#Sort the list in decreasing order 
genenrich = sort(genenrich, decreasing = TRUE)

#Now we can run the gene enrichment analysis
gse = gseGO(geneList=genenrich, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr",
             eps = 0)

##Visualization of enrichment
require(DOSE)
##dot plot
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
##ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

####################################################
########### Bonus: KEGG Pathway analysis ###########
####################################################

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids=bitr(names(genenrich), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
#Remove duplicate ids again
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

#Create a new vector which has only the genes which were successfully mapped using the bitr function above
kegg_gene_list = genenrich[names(genenrich) %in% dedup_ids$ENSEMBL] 

#Name vector with ENTREZ ids
names(kegg_gene_list) = dedup_ids$ENTREZID

#Omit any NA values 
kegg_gene_list=na.omit(kegg_gene_list)

#Sort the list in decreasing order again
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

##Humans organism again
kegg_organism = "hsa"

##kegg pathway analysis
kk2 = gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid",
               eps = 0)


#Produce the KEGG plott (look for the png file written by the below function on your computer)
#Choose "pathway.id" to the id of your interest. Look at kk2 to extract relevant ones.
dme = pathview(gene.data=kegg_gene_list, pathway.id=..., species = kegg_organism)
