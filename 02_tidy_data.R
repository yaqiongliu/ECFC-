

rm(list = ls())

load("output/exprSet_vst.Rdata")
load("metadata.Rdata")
load("DEseq2_DM_Diff.Rdata")

# 

colnames(metadata) <- c("sample","group")

### 
library(dplyr)
ensemble_symbol = res %>% 
  dplyr::select(gene_id,gene) %>% 
  filter(gene !="") 
rownames(exprSet_vst)<-str_replace_all(rownames(exprSet_vst),pattern = "\\.\\d+",replacement = "")

exprSet <- cbind(gene_id=rownames(exprSet_vst),exprSet_vst)

### 
exprSet <- merge(ensemble_symbol,exprSet,by="gene_id")

###

exprSet <- exprSet %>% 
  dplyr::select(-gene_id) %>% 
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(gene,.keep_all = T) %>% 
  dplyr::select(-newcolumn)

### 
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
### 
save(exprSet,file = "output/exprSet_DM1_heatmap.Rdata")

### 
exprSet <- t(exprSet)
exprSet <- as.data.frame(exprSet)


###
exprSet <- cbind(group= metadata$group,exprSet)

save(exprSet,file = "output/exprSet_DM1_tidy.Rdata")

###########################################################

rm(list = ls())
load(file = "output/exprSet_DM1_tidy.Rdata")

## steal plot
my_comparisons <- list(
  c("Diabetic", "Control")
)
library(ggpubr)
ggboxplot(
  exprSet, x = "group", y = "CXCR4",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")
library(dplyr)
test<-read.csv(file="different_genes.csv",header=T)
test<-subset(test,adj.P.Val<0.05)
test<-na.omit(test)
test<-test[order(test$logFC,decreasing = TRUE),]

