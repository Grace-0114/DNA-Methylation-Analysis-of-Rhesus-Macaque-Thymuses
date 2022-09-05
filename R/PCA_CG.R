system('g++ -v')
rm(list=ls())
library(ggbiplot)
library(plyr)
library(scales)
library(grid)
library(ggplot2)

setwd("D:/Rdata")
data = read.table(file ="D:/Rdata/DNA_methylation/result2200NAX.txt",head = T,sep ="\t",row.names = 1)
head(data)
data=as.data.frame(data,row.names = 1)
pca<-prcomp(t(data))


p<-ggbiplot(pca,var.axes = F,
         obs.scale = 1,
         ellipse = F,
         group=c('AA_1','AA_1','AA_1','AA_1','AB_1','AB_1','AB_1','AB_1','YA_9','YA_9','YA_9','YA_9','YB_9','YB_9','YB_9','YB_9'),
         labels =c('AA_9.gene','AA_9.promotor','AA_9.intron','AA_9.exon',
                   'AB_9.gene','AB_9.promotor','AB_9.intron','AB_9.exon',
                   'YA_1.gene','YA_1.promotor','YA_1.intron','YA_1.exon',
                   'YB_1.gene','YB_1.promotor','YB_1.intron','YB_1.exon'),labels.size = 3,scale = 1)
p
p + theme(panel.background = element_rect(fill = 'white', colour = 'black'))


p<-ggbiplot(pca,var.axes = F,
            obs.scale = 1,
            ellipse = F,
            group=c('YA_1','YA_1','YA_1','YA_1','YB_1','YB_1','YB_1','YB_1','AA_9','AA_9','AA_9','AA_9','AB_9','AB_9','AB_9','AB_9'),
            scale = 1)

setwd("D:/Rdata")
data = read.table(file ="D:/Rdata/DNA_methylation/result2199NAX.txt",head = T,sep ="\t",row.names = 1)
head(data)
data=as.data.frame(t(data),row.names = 1)
pca<-prcomp(data)


p<-ggbiplot(pca,var.axes = F,
            obs.scale = 1,
            ellipse = F,
            group=c('AA_1','AA_1','AA_1','AA_1','AB_1','AB_1','AB_1','AB_1','YA_9','YA_9','YA_9','YA_9','YB_9','YB_9','YB_9','YB_9'),
            labels =c('AA_9.gene','AA_9.promotor','AA_9.intron','AA_9.exon',
                      'AB_9.gene','AB_9.promotor','AB_9.intron','AB_9.exon',
                      'YA_1.gene','YA_1.promotor','YA_1.intron','YA_1.exon',
                      'YB_1.gene','YB_1.promotor','YB_1.intron','YB_1.exon'),labels.size = 3,scale = 1)
p
p + theme(panel.background = element_rect(fill = 'white', colour = 'black'))+coord_cartesian(xlim = c(-10, 25), ylim = c(-10, 10))
  
  
  scale_y_continuous(name = "This is y axis", breaks = seq(2, 4, by = 1), labels = c("-4", "-2", "0", "2", "4"))
