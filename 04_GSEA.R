
rm(list = ls())

##############################################################
### GSEA 
load(file = "DEseq2_ESR1_Diff.Rdata")
library(dplyr)
gene_df <- res %>% 
  dplyr::select(gene_id,logFC,gene) %>% 
  ## 
  filter(gene!="") %>% 
  ## 
  distinct(gene,.keep_all = T)

### 
geneList <- gene_df$logFC
### 
names(geneList) = gene_df$gene
## 
geneList = sort(geneList, decreasing = TRUE)

head(geneList)
library(clusterProfiler)

##
hallmarks <- read.gmt("resource/h.all.v7.1.symbols.gmt")
### 
gseahallmarks <- GSEA(geneList,TERM2GENE =hallmarks)

library(ggplot2)
yd <- as.data.frame(gseahallmarks)

dotplot(gseahallmarks,showCategory=30,split=".sign")+facet_grid(~.sign)

library(enrichplot)
pathway.id = "HALLMARK_ESTROGEN_RESPONSE_LATE"
gseaplot2(gseahallmarks, 
          color = "red",
          geneSetID = pathway.id,
          pvalue_table = T)

pathway.id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
gseaplot2(gseahallmarks, 
          color = "red",
          geneSetID = pathway.id,
          pvalue_table = T)

