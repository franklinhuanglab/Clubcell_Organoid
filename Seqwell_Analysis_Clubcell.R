

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggpubr)
library(clustree)
library(celda)
library(SingleCellExperiment)
theme_set(theme_cowplot())


#####
dge<-readRDS("./Organoid_decontX.rds")

####Separate the cells based on AR expression
dge_tmp<-subset(dge,subset=AR>0)
dge$AR_status="AR_Neg"
dge$AR_status[colnames(dge) %in% colnames(dge_tmp)]="AR_Pos"
dge$AR_ID<-paste0(dge$AR_status,"_",dge$ID)
dge<-SetIdent(dge,value = as.vector(dge$AR_ID))
dge@active.ident<- factor(dge@active.ident, levels = c("AR_Pos_Club","AR_Neg_Club",
                                                       "AR_Pos_BE","AR_Neg_BE",
                                                       "AR_Pos_Hillock","AR_Neg_Hillock",
                                                       "AR_Pos_Dividing","AR_Neg_Dividing"))

####Separate the cells based on treatment condition
dge$condition<-"NODHT"
dge$condition[dge$orig.ident %in% c("HNW_PR5269_T_org_DHT_V1_1_S26", "HNW_PR5269_T_org_DHT_V1_2_S27",
                                    "HNW_PR5316_T_apex_DHT_V1_1_S32", "HNW_PR5316_T_apex_DHT_V1_2_S33")]="DHT"
dge$condition[dge$orig.ident %in% c("HNW_PR5269_T_org_Enza_V1_1_S30", "HNW_PR5269_T_org_Enza_V1_2_S31",
                                    "HNW_PR5316_T_apex_Enza_V1_1_S36", "HNW_PR5316_T_apex_Enza_V1_2_S37")]="Enza"

dge$ID_Condition<-paste0(dge$ID,"_",dge$condition)


####Evaluate Module Score
Featurename<-c("BE","Hillock","Club","LE")
Features<-list()
BE_markers <- as.vector(read.csv("./Signature/BE_Markers.txt",header=F,sep="\n")$V1)
Hillock_markers <- as.vector(read.csv("./Signature/Hillock.txt",header=F,sep="\n")$V1)
Club_markers <- as.vector(read.csv("./Signature/Club_Markers.txt",header=F,sep="\n")$V1)
LE_markers <- as.vector(read.csv("./Signature/LE_Markers.txt",header=F,sep="\n")$V1)

###Make heatmap of canonical markers or customized genesets
library(tidyr)
library(textshape)
genelist<-c(BE_markers[2:21],Club_markers[2:21],Hillock_markers[2:21],ERGpos_markers,ERGneg_markers)
genelist<-unique(top10$gene)
dge <- ScaleData(object=dge,features=rownames(dge))
dge<-SetIdent(dge,value = as.vector(dge$condition))
dge_temp<-subset(dge,downsample=200)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$condition))

genelist<-c("PIGR","CP","MMP7","LTF","SCGB1A1","LCN2","PSCA","SCGB3A1","KLK3","KLK2","ACPP","NKX3-1","AR","KRT8","KRT18",
            "KRT5","KRT15","KRT17","KRT13")

library(tidyr)
z <- DoHeatmap(dge, features = genelist,assay = "RNA") + NoLegend()
y <- z$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()

my_sample_col <- data.frame((dge$condition),(dge$AR_status))
colnames(my_sample_col)=c("Condition","AR_Status")
my_gene_col <-as.data.frame(c(rep("Club_marker",8),rep("Luminal_marker",7),rep("Basal_Hillock_Tumor",4)))
rownames(my_gene_col)=t(make.unique(genelist))
colnames(my_gene_col)="Markers"
ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = F,clustering_method="ward.D2",cluster_cols = T,cluster_rows = T,scale = "none",show_rownames = T)
ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,annotation_row = my_gene_col,show_colnames = F,clustering_method="ward.D2",cluster_cols = T,cluster_rows = T,scale = "none",show_rownames = T)


####AR Signature analysis
dge<-SetIdent(dge,value = as.vector(dge$condition))
Hallmark_AR <- as.vector(read.csv("./Signature/HALLMARK_ANDROGEN_RESPONSE.txt",header=F,sep="\n")$V1)
Nelson_AR <- as.vector(read.csv("./Signature/NELSON_RESPONSE_TO_ANDROGEN_UP.txt",header=F,sep="\n")$V1)

dge<-AddModuleScore(object = dge,features = list(Hallmark_AR), ctrl = 5,name = "Hallmark_AR")
dge<-AddModuleScore(object = dge,features = list(Nelson_AR), ctrl = 5,name = "Nelson_AR")
dge@active.ident<- factor(dge@active.ident, levels = c("NODHT","DHT","Enza"))

dge$detailed_status<-paste0(dge$AR_status,"_",dge$condition)
dge<-SetIdent(dge,value = as.vector(dge$detailed_status))
dge@active.ident<- factor(dge@active.ident, levels = c("AR_POS_NODHT","AR_POS_DHT","AR_POS_Enza","AR_NEG_NODHT","AR_NEG_DHT","AR_NEG_Enza"))
Featurename="Nelson_AR"
V1<-as.vector(FetchData(object = dge,vars = paste0(Featurename,"1")))
BOX_df<-NULL
BOX_df$id<-dge@active.ident
BOX_df$value<-as.vector(V1[,1])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_violin() + geom_boxplot()  +#geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
  geom_signif(comparisons = list(c("AR_POS_NODHT","AR_POS_DHT")),map_signif_level=TRUE,y_position = 0.5, test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_POS_NODHT","AR_POS_Enza")),map_signif_level=TRUE,y_position = 0.55, test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_NEG_NODHT","AR_NEG_DHT")),map_signif_level=TRUE,y_position = 0.5, test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_NEG_NODHT","AR_NEG_Enza")),map_signif_level=TRUE,y_position = 0.55, test = "wilcox.test") + 
  #geom_signif(comparisons = list(c("NODHT","Enza")),map_signif_level=TRUE,y_position = 0.5, test = "wilcox.test") + 
  #geom_signif(comparisons = list(c("DHT","Enza")),map_signif_level=TRUE,y_position = 0.6, test = "wilcox.test") + 
  ggtitle(Featurename)
ggsave(file=paste0("Vlnplot_",Featurename,".pdf"),width = 30,height = 10,units = "cm")

#
###AR status
dge<-SetIdent(dge,value = as.vector(dge$AR_status))
Nelson_AR <- as.vector(read.csv("/Users/hsong/Desktop/Seqwell_combined/combined_analysis/geneset/NELSON_RESPONSE_TO_ANDROGEN_UP.txt",header=F,sep="\n")$V1)
Nelson_AR<-intersect(Nelson_AR,rownames(dge))
for (i in 1:length(Nelson_AR)){
V1<-as.vector(FetchData(object = dge,vars = Nelson_AR[i]))
BOX_df<-NULL
BOX_df$id<-dge@active.ident
BOX_df$value<-as.vector(V1[[1]])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(outlier.size = 0.5)  +#geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=4,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
  #geom_signif(comparisons = list(c("AR_POS","AR_NEG")),map_signif_level=TRUE,y_position = max(BOX_df$value), test = "wilcox.test") + NoLegend()+ 
  geom_signif(comparisons = list(c("AR_Pos_Club","AR_Neg_Club")),map_signif_level=TRUE,y_position = max(BOX_df$value), test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_Pos_BE","AR_Neg_BE")),map_signif_level=TRUE,y_position = max(BOX_df$value), test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_Pos_Hillock","AR_Neg_Hillock")),map_signif_level=TRUE,y_position = max(BOX_df$value), test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_Pos_Dividing","AR_Neg_Dividing")),map_signif_level=TRUE,y_position = max(BOX_df$value), test = "wilcox.test") + 
  ggtitle(Nelson_AR[i])+rotate_x_text(angle = 45)
ggsave(file=paste0("Vlnplot_",Nelson_AR[i],".pdf"),width = 10,height = 5)
}




Featurename="Hallmark_AR"
V1<-as.vector(FetchData(object = dge,vars = paste0(Featurename,"1")))
BOX_df<-NULL
BOX_df$id<-dge@active.ident
BOX_df$value<-as.numeric(V1[[1]])
BOX_df<-data.frame(BOX_df)
ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot()  +#geom_violin() + #geom_jitter(shape=16,position = position_jitter(0.1))+
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
  geom_signif(comparisons = list(c("AR_Pos_Club","AR_Neg_Club")),map_signif_level=TRUE,y_position = 0.5, test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_Pos_BE","AR_Neg_BE")),map_signif_level=TRUE,y_position = 0.55, test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_Pos_Hillock","AR_Neg_Hillock")),map_signif_level=TRUE,y_position = 0.5, test = "wilcox.test") + 
  geom_signif(comparisons = list(c("AR_Pos_Dividing","AR_Neg_Dividing")),map_signif_level=TRUE,y_position = 0.55, test = "wilcox.test") + 
  #geom_signif(comparisons = list(c("NODHT","Enza")),map_signif_level=TRUE,y_position = 0.5, test = "wilcox.test") + 
  #geom_signif(comparisons = list(c("DHT","Enza")),map_signif_level=TRUE,y_position = 0.6, test = "wilcox.test") + 
  ggtitle(Featurename) +rotate_x_text(angle = 45)
ggsave(file=paste0("Vlnplot_",Featurename,".pdf"),width = 30,height = 15,units = "cm")




