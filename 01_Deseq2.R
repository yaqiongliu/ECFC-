rm(list = ls())

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(vsn)
library(pheatmap)
library(viridis)
library(PoiClaClu)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggrepel)
library(pcaExplorer)
library(msigdbr)
library(fgsea)
library(patchwork)
library(pathview)
library(KEGGprofile)
library(scales)
library(ggrepel)
load("exprSet.Rdata")
load("metadata.Rdata")

colnames(metadata) <- c("sample","group")
#########################################################

library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                             colData=metadata, 
                             design=~group)
nrow(dds)
rownames(dds)

dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

#########################################################

vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, "group")


exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:10,1:7]

save(exprSet_vst,file = "output/exprSet_vst.Rdata")

#########################################################

### estimating size factors
### estimating dispersions
### gene-wise dispersion estimates
### mean-dispersion relationship
### final dispersion estimates
### fitting model and testing
dds <- DESeq(dds)

#########################################################

contrast=c("group", "Diabetic", "Control")

dd1 <- results(dds, contrast=contrast, alpha = 0.05)

DESeq2::plotMA(dd1, ylim=c(-5,5))

dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
DESeq2::plotMA(dd2, ylim=c(-5,5))


#########################################################

library(dplyr)
library(tibble)
library(tidyr)
res <- dd1 %>% as.data.frame()%>% 
rownames_to_column(var = "transcripts") %>%
  arrange(pvalue) %>%
  mutate(ensembl = str_replace_all(transcripts,pattern = "\\.\\d+",replacement = ""))
  
#save(res,file = "res_anno_demo.Rdata")
#########################################################

library(AnnotationDbi)
library(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=res$ensembl,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
                 res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$ensembl,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")    
                     
###
res<-res[,-1]
res<-res[,-4]

colnames(res) <- c("baseMean","logFC","lfcSE","P.Value","adj.P.Val","gene_id","gene","entrez")

#########################################################
save(res,file = "DEseq2_DM_Diff.Rdata")
