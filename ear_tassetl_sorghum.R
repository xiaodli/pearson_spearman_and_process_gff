### single anchorwave
library(ggplot2)
data =read.table("/home/xiaodong/Desktop/recent/correlation/anchorwave_R_corr.csv", 
                 header = TRUE, sep=",")
# new_df <- data[c("corr_value1", "corr_value2")]
data$tas_ear = factor(data$tas_ear, levels=c("tassel", "ear"))
plot <- ggplot(data, aes(x=tas_ear, y=corr_value, fill=tas_ear)) + geom_boxplot() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  labs(fill = "tassel_ear")
  # geom_jitter(color="blue", size=0.1, alpha=0.2)
  # facet_grid(~tas_ear)
png("/home/xiaodong/Desktop/recent/correlation/result/anchorwave.corr.png", 
    width=547, height=394, res=120)
plot
dev.off()


### anchorwave and wgdi
library(ggplot2)
data =read.table("/home/xiaodong/Desktop/recent/correlation/anchorwave_wgdi_corr_R.csv", 
                 header = TRUE, sep=",")
# new_df <- data[c("corr_value1", "corr_value2")]
data$le_ri = factor(data$le_ri, levels=c("AnchorWave", "WGDI"))
plot <- ggplot(data, aes(x=tas_ear, y=corr_value, fill=tas_ear)) + geom_boxplot() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + facet_grid(~le_ri) +
  labs(fill = "tassel_ear")
# geom_jitter(color="blue", size=0.1, alpha=0.2)

png("/home/xiaodong/Desktop/recent/correlation/result/anchorwave.wgdi.corr.png", 
    width=547, height=394, res=120)
plot
dev.off()


### wgdi diff
library(ggplot2)
data =read.table("/home/xiaodong/Desktop/recent/correlation/wgdi.corr.diff.R.csv", 
                 header = TRUE, sep=",")
# new_df <- data[c("corr_value1", "corr_value2")]
data$le_ri = factor(data$le_ri, levels=c("same", "diff"))
plot <- ggplot(data, aes(x=tas_ear, y=corr_value, fill=tas_ear)) + geom_boxplot() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + facet_grid(~le_ri) +
  labs(fill = "tassel_ear")
# geom_jitter(color="blue", size=0.1, alpha=0.2)

png("/home/xiaodong/Desktop/recent/correlation/result/wgdi.diff.corr.png", 
    width=547, height=394, res=120)
plot
dev.off()


### anchorwave and wgdi diff are merged
library(ggplot2)
data =read.table("/home/xiaodong/Desktop/recent/correlation/merge.csv", 
                 header = TRUE, sep=",")
# new_df <- data[c("corr_value1", "corr_value2")]
data$le_ri = factor(data$le_ri, levels=c("AnchorWave", "WGDI", "same", "diff"))
plot <- ggplot(data, aes(x=tas_ear, y=corr_value, fill=tas_ear)) + geom_boxplot() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + facet_grid(~le_ri) +
  labs(fill = "tassel_ear")
# geom_jitter(color="blue", size=0.1, alpha=0.2)

png("/home/xiaodong/Desktop/recent/correlation/result/merge.corr.png", 
    width=1047, height=394, res=120)
plot
dev.off()


### anchorwave and wgdi block merged
library(ggplot2)
# data =read.table("/home/xiaodong/Desktop/recent/correlation/anchorwave_wgdi_block.csv", 
#                  header = TRUE, sep=",")
data =read.table("/home/xiaodong/Desktop/recent/correlation/anchorwave_wgdi_block5555555555.csv",
                 header = TRUE, sep=",")
data$le_ri = factor(data$le_ri, levels=c("AnchorWave", "WGDI"))
plot <- ggplot(data, aes(x=tas_ear, y=corr_value, fill=tas_ear)) + geom_boxplot() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + facet_grid(~le_ri) +
  labs(fill = "tassel_ear")
png("/home/xiaodong/Desktop/recent/correlation/result/anchorwave_wgdi_block.corr5555.png", 
    width=594, height=394, res=120)
plot
dev.off()



## 
### left inner right
library(ggplot2)
data =read.table("/home/xiaodong/Desktop/recent/correlation/correlation_inner_left_right.csv", 
                 header = TRUE, sep=",")
# new_df <- data[c("corr_value1", "corr_value2")]
data$tas_ear = factor(data$tas_ear, levels=c("tassel", "ear"))
data$source = factor(data$source, levels=c("inner", "AnchorWave", "WGDI"))
plot <- ggplot(data, aes(x=source, y=corr_value, fill=source)) + geom_boxplot() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size = 7),
        legend.key.size = unit(0.5, "cm")) + facet_grid(~tas_ear) +
  labs(fill = "source")
# geom_jitter(color="blue", size=0.1, alpha=0.2)

png("/home/xiaodong/Desktop/recent/correlation/result/left.inner.right.corr.png", 
    width=547, height=394, res=120)
plot
dev.off()



## venn 
library(cowplot)
library(VennDiagram)
data1 =read.table("/home/xiaodong/Desktop/recent/correlation/result/correlation.csv", 
                 header = TRUE, sep=",")
data2 =read.table("/home/xiaodong/Desktop/recent/correlation/result/wgdi_correlation.csv", 
                 header = TRUE, sep=",")
Anchor_tassel <- unique(na.omit(data1$tassel_gene_name))
Anchor_ear <- unique(na.omit(data1$ear_gene_name))
WGDI_tassel <- unique(na.omit(data2$tassel_gene_name))
WGDI_ear <- unique(na.omit(data2$ear_gene_name))
tassel_list <- list(AnchorWave_tassel=Anchor_tassel, WGDI_tassel=WGDI_tassel)
ear_list <- list(AnchorWave_ear=Anchor_ear, WGDI_ear=WGDI_ear)

tassel <- venn.diagram(x=tassel_list, filename = NULL, alpha=c(0.6, 0.6), 
                       euler.d=T, scaled=T, fill=c("#0073C2FF",colors()[468]),
                       lwd=c(1,1),cex=1,
                       cat.pos=c(-24,12),
                       resolution = 500,
                       cat.dist=c(0.05,0.02),
                       cat.cex=1)
cowplot::plot_grid(tassel)
# png("/home/xiaodong/Desktop/tassel.png", 
    # width=800, height=459, res=500)
# tassel
# dev.off()




library(cowplot)
library(VennDiagram)
data1 =read.table("/home/xiaodong/Desktop/recent/correlation/result/correlation.csv", 
                  header = TRUE, sep=",")
data2 =read.table("/home/xiaodong/Desktop/recent/correlation/result/wgdi_correlation.csv", 
                  header = TRUE, sep=",")
Anchor_tassel <- unique(na.omit(data1$tassel_gene_name))
Anchor_ear <- unique(na.omit(data1$ear_gene_name))
WGDI_tassel <- unique(na.omit(data2$tassel_gene_name))
WGDI_ear <- unique(na.omit(data2$ear_gene_name))
tassel_list <- list(AnchorWave_tassel=Anchor_tassel, WGDI_tassel=WGDI_tassel)
ear_list <- list(AnchorWave_ear=Anchor_ear, WGDI_ear=WGDI_ear)
ear <- venn.diagram(x=ear_list, filename = NULL, alpha=c(0.6, 0.6), 
                       euler.d=T, scaled=T, fill=c("#0073C2FF",colors()[468]),
                       lwd=c(1,1),cex=1,
                       cat.pos=c(-24,12),
                       resolution = 500,
                       cat.dist=c(0.05,0.02),
                       cat.cex=1)
cowplot::plot_grid(ear)
# png("/home/xiaodong/Desktop/ear.png", width=800, height=459, res=500)
# ear
# dev.off()

