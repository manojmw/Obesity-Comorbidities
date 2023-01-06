#!/usr/bin/env Rscript

# importing library
library(clusterProfiler)
library(org.Hs.eg.db)

library("DOSE")
library("GO.db")
library("GSEABase")


library("dplyr")
library("tidyr")
library("ggplot2")
library("stringr")
library("RColorBrewer")
library("rWikiPathways")
library("enrichplot")
library("Rgraphviz")
library("topGO")


# Input data - Hub gene list
CRa <- read.table("hubgenes.txt", header=T)


#### For only enrichment
CRa_en <- CRa[,1]

#### GO over-representation analysis
CRa_ego <- enrichGO(gene      = CRa_en,
				keyType		  = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "none",
                readable      = TRUE)

#### GO over-representation analysis for MF
CRa_ego_MF <- enrichGO(gene   = CRa_en,
				keyType		  = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "none",
                readable      = TRUE)
### GO over-representation analysis for CC
CRa_ego_CC <- enrichGO(gene   = CRa_en,
				keyType		  = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "none",
                readable      = TRUE)
### GO over-representation analysis for BP
CRa_ego_BP <- enrichGO(gene   = CRa_en,
				keyType		  = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "none",
                readable      = TRUE)

# Generate plots and write to file
pdf("MF_dot.pdf")
dotplot(CRa_ego_MF, showCategory = 20)
dev.off()

pdf("CC_dot.pdf")
dotplot(CRa_ego_CC, showCategory = 20)
dev.off()

pdf("BP_dot.pdf")
dotplot(CRa_ego_BP, showCategory = 20)
dev.off()

pdf("BP_dot.pdf")
dotplot(CRa_ego, showCategory = 20)
dev.off()

write.table(CRa_ego,file="ALL_enrichment.txt", sep="\t")


####################################
####################################			 
### Tree plot
edox2 <- pairwise_termsim(CRa_ego)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')


####################################
####################################


#### KEGG enrichment 
CRa_gene_entrez = bitr(CRa_en, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
CRa_gene_entrez <- CRa_gene_entrez[,2]
CRa_kk <- enrichKEGG(gene      = CRa_gene_entrez,
                 organism      = 'hsa',
				 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)

# Generate plots and write to file	 
write.table(CRa_kk,file="KEGG_enrichment.txt", sep="\t")

pdf("KEGG_barplot.pdf")
barplot(CRa_kk, showCategory = 20)
dev.off()

#### Wikipathway
# enrichWP(CRa_gene_entrez_uf, organism = "Mus musculus") 


#### Over-representation analysis for disease ontology
x <- enrichDO(gene          = CRa_gene_entrez,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
			  

# Generate plots and write to file	  
pdf("Disease_barplot.pdf")
barplot(x, showCategory = 20)
dev.off()


#### For displaying only top 10 terms
mutate(CRa_ego, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore")
	
#### Two coulmn format

### Modify accoing to the format for GSEA
CRa_gsea <- CRa[,2]
names(CRa_gsea) <- as.character(CRa[,1])
CRa_gsea <- sort(CRa_gsea, decreasing = TRUE)

	
### Gene-Concept Network
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

#For enriched terms
p1 <- cnetplot(CRa_ego, foldChange=CRa_gsea)

## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(CRa_ego, categorySize="pvalue", foldChange=CRa_gsea)
p3 <- cnetplot(CRa_ego, foldChange=CRa_gsea, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

#For KEGG pathway terms 
CRa_kkx <- setReadable(CRa_kk, 'org.Hs.eg.db', 'ENTREZID')


p1 <- cnetplot(CRa_kkx, foldChange=CRa_gsea)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(CRa_kkx, categorySize="pvalue", foldChange=CRa_gsea)
p3 <- cnetplot(CRa_kkx, showCategory = 20, foldChange=CRa_gsea, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

write.table(CRa_kkx,file="KEGG_enrichment.txt", sep="\t")

# For Disease enriched terms 
CRa_dix <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')


p1 <- cnetplot(CRa_dix, foldChange=CRa_gsea)
## categorySize can be scaled by 'pvalue' or 'geneNum'

p2 <- cnetplot(CRa_dix, categorySize="pvalue", foldChange=CRa_gsea)
p3 <- cnetplot(CRa_dix, showCategory = 20, foldChange=CRa_gsea, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

#### Tree plot
# For enriched terms
CRa_ego2 <- pairwise_termsim(CRa_ego)
p1 <- treeplot(CRa_ego2)
p2 <- treeplot(CRa_ego2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

#### Heatplot
p1 <- heatplot(CRa_kkx, showCategory=20)