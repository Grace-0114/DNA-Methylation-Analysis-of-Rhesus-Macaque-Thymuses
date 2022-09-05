#2019-10-15#
#Grace#
#甲基化各样本基因count相关系数#

rm(list=ls())
library(pheatmap)
data=read.table("C:/Users/pc/Desktop/result2200NA.txt",head = TRUE, row.names = 1)

matrix=cor(data, method="spearman") 
head(matrix)


# Create the heatmap annotation
colnames(matrix) <- c("A_9-genebody","A_9-promotor","A_9-intron","A_9-exon","B_9-genebody","B_9-promotor","B_9-intron","B_9-exon","A_1-genebody","A_1-promotor","A_1-intron","A_1-exon","B_1-genebody","B_1-promotor","B_1-intron","B_1-exon")#paste("Test", 1:16, sep = ",") 
rownames(matrix) <- c("A_9-genebody","A_9-promotor","A_9-intron","A_9-exon","B_9-genebody","B_9-promotor","B_9-intron","B_9-exon","A_1-genebody","A_1-promotor","A_1-intron","A_1-exon","B_1-genebody","B_1-promotor","B_1-intron","B_1-exon")#"genebody", 1:16, sep = ",")  

annotation_col = data.frame(Age = c("9","9","9","9","9","9","9","9","1","1","1","1","1","1","1","1"), 
                            Gender =c("F","F","F","F","M","M","M","M","F","F","F","F","M","M","M","M"))
rownames(annotation_col) = c("A_9-genebody","A_9-promotor","A_9-intron","A_9-exon","B_9-genebody","B_9-promotor","B_9-intron","B_9-exon","A_1-genebody","A_1-promotor","A_1-intron","A_1-exon","B_1-genebody","B_1-promotor","B_1-intron","B_1-exon")

annotation_row = data.frame(Region = c("genebody","promotor","intron","exon","genebody","promotor","intron","exon","genebody","promotor","intron","exon","genebody","promotor","intron","exon"))
rownames(annotation_row) = c("A_9.genebody","A_9.promotor","A_9.intron","A_9.exon","B_9.genebody","B_9.promotor","B_9.intron","B_9.exon","A_1.genebody","A_1.promotor","A_1.intron","A_1.exon","B_1.genebody","B_1.promotor","B_1.intron","B.1-exon")

ann_colors =list(Age = c("1"="#96D4B4","9"="#409B6C"),
                 Gender = c(F="#E47AA1",M="#6CA3B4"),
                 Region = c(genebody="#FAE3CA",promotor="#F2C9BA",intron="#D1A380",exon="#A7C2C2"))

pheatmap(matrix,cluster_rows = TRUE,
         cluster_cols = TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         cellwidth = 15, cellheight = 12,angle_col = 45,
         annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = ann_colors,filename = "M_cor.pdf")