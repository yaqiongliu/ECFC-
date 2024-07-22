
rm(list = ls())

load(file = "output/exprSet_heatmap.Rdata")

load(file = "DEseq2_DM1_Diff.Rdata")
### 
library(dplyr)
diffgene <- res %>% 
  filter(gene !="") %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1)
save(diffgene,file = "output/diffgene.Rdata")


library(pheatmap)
### 
heatdata <- exprSet[diffgene$gene,]


load(file = "metadata.Rdata")

annotation_col <- data.frame(group=metadata$group)
rownames(annotation_col) <- metadata$sample

### 
pheatmap(heatdata)

### 
pheatmap(heatdata, 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col =annotation_col, 
         annotation_legend=TRUE,  
         show_rownames = F,
         show_colnames = F,
         scale = "row", 
         color =colorRampPalette(c("blue", "white","red"))(100),
         #filename = "heatmap_F.pdf",
         cellwidth = 25, cellheight = 0.5,
         fontsize = 10)


################################################################

library(ggplot2)
library(ggrepel)
data <- res
data$significant <- as.factor(data$P.Value<0.05 & abs(data$logFC) > 0.5)
ggplot(data=data, aes(x=logFC, y =-log10(P.Value),color=significant)) +
  geom_point(alpha=0.8, size=1.2,col="black")+
  geom_point(data=subset(data, logFC > 2),alpha=0.8, size=1.4,col="red")+
  geom_point(data=subset(data, logFC < -2),alpha=0.6, size=1.4,col="blue")+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(-2,2),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  geom_point(data=subset(data, abs(logFC) >= 4),alpha=0.8, size=3,col="green")+
  geom_text_repel(data=subset(data, abs(logFC) > 4), 
                  aes(label=gene),col="black",alpha = 0.8)
