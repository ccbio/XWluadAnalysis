library(ggvenn)
library(igraph)
library(tidyr)# 使用的gather & spread
library(reshape2) # 使用的函数 melt & dcast 
library(dplyr)#mutate

pdata_LUAD=read.csv(sep='\t','../../../../download/xenabrowser/LUAD/GDC TCGA Lung Adenocarcinoma (LUAD)/TCGA-LUAD.GDC_phenotype.tsv')
pdata_LUSC=read.csv(sep='\t','../../../../download/xenabrowser/LUSC/GDC TCGA Lung Squamous Cell Carcinoma (LUSC)/TCGA-LUSC.GDC_phenotype.tsv')
pdata_EAS=read.csv(sep='\t','../../zww东亚/e2211c4ba45274f8a4763217c0b0d20e/GIS031_GSK/GIS031.clinical.patient.txt')
tumor_XW=read.csv('../../group.csv')
tumor_XW$name=paste0(tumor_XW$batch,'_',gsub('-','\\.',tumor_XW$filename))
##构建表达矩阵
pdata_LUAD=pdata_LUAD[,c("submitter_id.samples",
"sample_type.samples",
'gender.demographic','cigarettes_per_day.exposures')]
colnames(pdata_LUAD)[1]="Tumor_Sample_Barcode"
pdata_LUSC$smoker='YES'
pdata_LUAD$smoker[is.na(pdata_LUAD$cigarettes_per_day.exposures)]='NO'
NoFemale_TCGA=pdata_LUAD[pdata_LUAD$gender.demographic=='female'&pdata_LUAD$smoker=='NO',]
NoFemale_TCGA=merge(NoFemale_TCGA,pdata_LUAD)

NoFemale_EAS=pdata_EAS[pdata_EAS$GENDER=='Female'&pdata_EAS$SMOKING_STATUS=='No',]



target_sample=intersect(NoFemale_TCGA$Tumor_Sample_Barcode,colnames(tcga_lusc2))
target_NoFemale_TCGA=NoFemale_TCGA[NoFemale_TCGA$Tumor_Sample_Barcode %in% target_sample,]
matrix=tcga_lusc4[,target_sample]
{
  x=gene_trans[gene_trans$Gene.type=='protein-coding',]
  x=na.omit(x)
  XW_luad3=XW_luad2[x$Ensembl.ID,]
  XW_luad3$gene_id=x[,2]
  XW_luad4=aggregate(.~gene_id,data=XW_luad3,sum)
  row.names(XW_luad4)=XW_luad4[,1]
  XW_luad4=XW_luad4[,-1]
}
colnames(EAS_luad_normal2)=paste0(colnames(EAS_luad_normal2),'_normal')
EAS_luad_all=cbind(EAS_luad_normal2,EAS_luad2)

os=c(rep(-1,88),rep(1,172))



meta=tumor_XW#样本信息
row.names(meta)=meta[,1]
meta$typel=NA
meta$typel[meta$GROUP=='normal']<- -1
meta$typel[meta$GROUP=='tumor']<- 1
matrix=as.data.frame(XW_luad2)
#NJ=gd1_wide2#NJ
matrix=t(XW_luad4)#按照NJ的基因顺序调整表达矩阵，列为NJ的基因顺序,行为样本
matrix=cbind(meta$typel,matrix[intersect(meta$name,row.names(matrix)),])

V1=c(rep(-1,88),rep(1,172))
matrix=cbind(V1,t(EAS_luad_all4))
saveRDS(tcga,file="../WES-1/STSVM/matrix_XW.rdata")


##Y列表（normal vs. tumor）
matrix=readRDS("../../WES/matrix_LUAD.rdata")
matrix=readRDS("matrix_XW.rdata")
matrix=readRDS("matrix_EAS.rdata")
#stsvm
##构建邻接矩阵
matrix=as.data.frame(t(tcga[,-1]))
ppi=read.table('../../WES/STSVM/PathwayCommons12.All.hgnc.sif')
gene=as.data.frame(rownames(matrix))
pp1=ppi
colnames(pp1)=c("gene","cor","gene2")
colnames(gene)="gene"
pp2=merge(pp1,gene,by=1)
ppi1=pp2[which(pp2$cor=='interacts-with'),]
ppi2=pp2[which(pp2$cor=='controls-state-change-of'),]
ppi3=rbind(ppi1,ppi2)
ppi3$var=1
ppi4=ppi3[,-2]
library(igraph)
library(tidyr)# 使用的gather & spread
library(reshape2) # 使用的函数 melt & dcast 
library(dplyr)#mutate
mygraph <- graph.data.frame(ppi4[,c(1,2)])
gd1_wide1=get.adjacency(mygraph, sparse = FALSE)
saveRDS(gd1_wide1,file = "XW.adj.matrix.rdata")

##构建表达矩阵和Y列表（y取xw和
#rm(gd1_wide1)
matrix=readRDS("../../WES/matrix_LUAD.rdata")
matrix=readRDS("matrix_XW.rdata")
matrix=readRDS("matrix_EAS.rdata")
matrix=as.data.frame(matrix)
matrix=matrix[order(matrix$V1),]
##构建re列表
re=matrix[,1]
re=as.data.frame(re)
re$re=as.factor(re$re)
re=as.list(re)
##构建express矩阵
exp=matrix[,-1]
exp=as.data.frame(exp)
exp[,1:18778]=lapply(exp[,1:18778],as.numeric)

f<-function(x) sum(x==0)
#apply(exp,2,f)
n0 <- apply(exp,2,f)
# 统计0元素大于2个的行号
i0 <- which(n0 > 150*0.9)
# 删除0元素大于2个的行
exp=exp[,-i0 ]




exp1=as.matrix.data.frame(exp)
exp1=list(exp1)
list=c(exp1,re)
names(list)=c("genes","y")
saveRDS(list,file="XW_LUAD.list.rdata")

##cvstsvm训练
library(netClass)
expr1=readRDS("XW_LUAD.list.rdata")
ad.matrix1=readRDS("XW.adj.matrix.rdata")
x=expr1$genes
y=expr1$y
r.stsvm1 <- cv.stsvm(x=x,x.mi=NULL,y=y,folds=10,Gsub=ad.matrix1,op.method="pt",
                     repeats=10, parallel=FALSE, cores=2,DEBUG=TRUE,
                     pt.pvalue=0.01,op=0.80,
                     aa=5,a=1,p=2,allF=TRUE, seed=1234,Cs=10^(-3:3))

saveRDS(r.stsvm1,'XW_LUAD.stsvm.RDS')
mean(r.stsvm1$auc)

# STSVM GENE-----------------------

genelist.stsvm=function(r.stsvm){
  genesstsvm=c()
  
  for (i in 1:10){
    repeatl=paste0("Repeat",i)
    for (i in 1:10){
      foldl=paste0("Fold",i)
      x=r.stsvm[["fits"]][[repeatl]][[foldl]][["trained"]][["features"]]
      genesstsvm=c(genesstsvm,x)
    }
  }
  return(genesstsvm)
}
xx=as.data.frame(LU_TCGA.stsvm$auc)
max(xx)

gene_res=function(df){
  LU_TCGA.stsvm=readRDS(df)
  LU_TCGA.genelist=as.data.frame(table(genelist.stsvm(LU_TCGA.stsvm)))
  LU_TCGA.gene_Res=as.data.frame(subset(LU_TCGA.genelist,LU_TCGA.genelist$Freq>=100))
 return(LU_TCGA.gene_Res)
}

LUAD_stsvm_gene=gene_res("TCGA_LUAD.stsvm.RDS")
EAS_stsvm_gene=gene_res("EAS_LUAD.stsvm.RDS")
XW_stsvm_gene=gene_res("XW_LUAD.stsvm.RDS")

GENE_list=list(TCGA=as.character(LUAD_stsvm_gene$Var1),
               EAS=as.character(EAS_stsvm_gene$Var1),
               XW=as.character(XW_stsvm_gene$Var1))
saveRDS(GENE_list,file = 'GENE_list.stsvm.rds')

ggvenn(
  data =GENE_list,         # 数据列表
  columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选
  show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
  label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
  show_percentage = T,      # 显示每一组的百分比
  digits = 1,               # 百分比的小数点位数
  fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"), # 填充颜色
  fill_alpha = 0.5,         # 填充透明度
  stroke_color = "white",   # 边缘颜色
  stroke_alpha = 0.5,       # 边缘透明度
  stroke_size = 0.5,        # 边缘粗细
  stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
  set_name_color = "black", # 组名颜色
  set_name_size = 6,        # 组名大小
  text_color = "black",     # 交集个数颜色
  text_size = 4             # 交集个数文字大小
)
# GO enrichment ---------------------
library(R.utils)
library(openxlsx)
library(ggplot2)#
library(stringr)
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)
library(ggvenn)
raw_beautifuk_dotplot <- function(raw_plot){
  data <- raw_plot$data
  # data$Description <- factor(data$Description,
  #                            levels = unique(data$Description))
  # # data$Cluster <- factor(data$Cluster, levels = rev(unique(data$Cluster)))
  p1 = ggplot(data,aes(Cluster,Description,size = GeneRatio))+
    geom_point(shape=21,aes(fill= qvalue),position =position_dodge(0))+
    theme_minimal()+xlab(NULL)+ylab(NULL) +
    scale_size_continuous(range=c(1,10))+
    theme_bw()+
    scale_fill_gradient(low = "#E54924", high = "#498EA4")+
    theme(legend.position = "right",legend.box = "vertical",
          legend.margin=margin(t= 0, unit='cm'),
          legend.spacing = unit(0,"in"),
          axis.text.x  = element_text(color="black",size=16),
          axis.text.y  = element_text(color="black",size=15),
          legend.text = element_text(size =12,color="black"),
          legend.title = element_text(size =12,color="black"))
  #   ) +
  #   scale_y_discrete(labels= data$Description)
  return(p1)
}
aging_cell_function_plot=function(a,b,c,DD){
  #a=unique(unlist(strsplit(paste(a$genes, collapse = ","),',')))
  #b=unique(unlist(strsplit(paste(b$genes, collapse = ","),',')))
  
  d=as.data.frame(cbind(a,b,c))
  colnames(d)=DD
  y = compareCluster(d, fun='enrichGO', 
                     OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont="ALL",
                     qvalueCutoff =0.05,pvalueCutoff =0.05)
  R.utils::setOption("clusterProfiler.download.method",'auto')
  e=d
  
  a <- bitr(a,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
  b <- bitr(b,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
  c <- bitr(c,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
  d=as.data.frame(cbind(a$ENTREZID,b$ENTREZID,c$ENTREZID))
  colnames(d)=DD
  y = compareCluster(d, fun='enrichKEGG', 
                     organism="hsa", #keyType = 'SYMBOL',
                     qvalueCutoff =0.01,pvalueCutoff =0.01)
  d=as.data.frame(cbind(b$ENTREZID,c$ENTREZID))
  colnames(d)=DD[c(2,3)]
  y = compareCluster(d, fun='enrichKEGG', 
                     organism="hsa", #keyType = 'SYMBOL',
                     qvalueCutoff =0.01,pvalueCutoff =0.01)
  KEGG<-enrichKEGG(gene=c$ENTREZID,
                   organism = "hsa",#keyType = "SYMBOL",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01)
  dotplot(KEGG,showCategory = 10)
  GO_plot <- dotplot(y, showCategory =10) + 
    theme(axis.text.x = element_text(angle=45, hjust=1))
  p=raw_beautifuk_dotplot(GO_plot)
  return(p)
}
p_XW_EAS=aging_cell_function_plot(a=GENE_list$XW,b=GENE_list$EAS,c=GENE_list$TCGA,DD=c('XW','EAS','TCGA'))
p_XW_LUAD=aging_cell_function_plot(a=XW_only_LUAD,b=LUAD_only_XW,DD=c('XW','TCGA'))

####################
matrix=readRDS("matrix_LUAD.rdata")
matrix=readRDS("matrix_XW.rdata")
matrix=readRDS("matrix_EAS.rdata")
matrix=as.data.frame(matrix)
matrix=matrix[order(matrix$V1),]
plot_matrix=matrix[,intersect(XW_stsvm_gene$Var1,colnames(matrix))]
#plot_matrix=matrix[,LUAD_stsvm_gene$Var1]

annotation_row = data.frame(
  type = as.character(matrix[,1])
)

row.names(annotation_row)=row.names(matrix)
####################microbiome
annotation_row = data.frame(
  Insomnia = as.character(grouppp1$sleep),
  Periodontitis=as.character(grouppp1$tooth)
)
row.names(annotation_row)=paste0('sub',grouppp1$ID)

annotation_col = data.frame(Pheno=gene_Res,row.names = gene_Res)
annotation_col$Pheno='insomnia'
annotation_col$Pheno[row.names(annotation_col) %in% tooth_metab]='periodontitis'
annotation_col$Pheno[row.names(annotation_col) %in% Both]='Both'
#row.names(annotation_col) <- annotation_col[,1]
#annotation_col=annotation_col[,2]
toothcolor <- c("red","#016D06") 
names(toothcolor) <- c('1','-1') #类型颜色

sleepcolor <- c('#F3B1A0','#68A180')
names(sleepcolor) <- c('1','-1')

targetcolor <- c("green","blue","orange")
names(targetcolor) <- c("Both","insomnia","periodontitis" )
ann_colors <- list(Periodontitis=toothcolor, Insomnia= sleepcolor,
                   Pheno=targetcolor) #颜色设置
matrix=readRDS("matrix_XW.rdata")

plot_matrix=na.omit(plot_matrix)
pheatmap::pheatmap(t(plot_matrix),
                   scale = 'row',
                   clustering_distance_cols = "correlation", 
                   #clustering_distance_cols = "euclidean",
                   clustering_method="mcquitty",
                   cluster_cols = T,
                   cluster_rows = F,
                   #gaps_col = c(6,12),
                   #gaps_row = c(6),
                   show_rownames = F,
                   show_colnames = F,
                   border_color = "lightgrey",
                   #annotation_row = annotation_row,
                   annotation_col = annotation_row#,
                   #annotation_colors = ann_colors
)

pca <- prcomp(plot_matrix, scale=T)
library(factoextra)
#pdf("1.pdf")
fviz_eig(pca, addlabels = TRUE)
#dev.off()
fviz_pca_ind(pca, repel=T) 
fviz_pca_ind(pca, col.ind=annotation_row$type, mean.point=F, 
             addEllipses = T, legend.title="Groups")

library(umap)
library(ggforce)
pca.umap=umap(pca$x)
plot_umap=pca.umap$layout
plot_umap1=cbind(matrix$V1,plot_umap)
colnames(plot_umap1)=c('type','UMAP_1','UMAP_2')

ggplot(plot_umap1)+
  geom_point(aes(x=UMAP_1,y=UMAP_2,color=as.factor(type)),size=5)+
  #stat_ellipse(aes(x=UMAP_1,y=UMAP_2,fill=as.factor(type)),alpha=1)+
  theme_bw()
  
#---画柱图带标准差 ---------------------------
library(tidyverse)
library(rstatix)
library(ggpubr)


test=cbind(grouppp2,log(stsvm_gene2+1))
test=test[,c('ID','group',colnames(stsvm_gene2))]
test$group=gsub("[0-9]","",test$group)
######metab_select
test=cbind(grouppp2,stsvm_gene)
test=test[,c('ID','group',colnames(stsvm_gene))]
test$group=gsub("[0-9]","",test$group)

set.seed(123)

test <- test %>%
  #gather(key = "time", value = "score", t1, t2, t3) %>%
  gather(key = "mircobiome_sub", value = "value",
         -group,-ID)%>%convert_as_factor(ID)

plot_micro<-function(M_name){
  #data=test1
  data=subset(test,test$mircobiome_sub==M_name)
  ggplot(data,aes(x=group,y=value))+#指定数据
    stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
    geom_boxplot(aes(fill=group), #绘制箱线图函数
                 outlier.colour="white",size=0.8)+#异常点去除
    theme(panel.background =element_blank(), #背景
          axis.line=element_line(),#坐标轴的线设为显示
          legend.position="none",plot.title = element_text(size=14))+#图例位置
    geom_jitter(width = 0.2)+#添加抖动点
    ggtitle(gsub('s__','',M_name))+#标题
    geom_signif(comparisons = list(c("CTR","PD"),
                                   c("PD","PS"),
                                   c("CTR","PS")),# 设置需要比较的组
                map_signif_level = T, #是否使用*显示
                test = t.test, ##计算方法
                y_position = c(10,11,12),#图中横线位置设置
                tip_length = c(c(0.01,0.01),
                               c(0.01,0.01),
                               c(0.01,0.01)),#横线下方的竖线设置
                size=0.8,color="black")+
    theme_prism(palette = "candy_bright",
                base_fontface = "plain", # 字体样式，可选 bold, plain, italic
                base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
                base_size = 12,  # 图形的字体大小
                base_line_size = 0.8, # 坐标轴的粗细
                axis_text_angle = 0)+ # 可选值有 0，45，90，270
    scale_fill_prism(palette = "candy_bright")
  #  return(p)
}

umap_list <- lapply(unique(test$mircobiome_sub), plot_micro)
plot_grid(plotlist = umap_list, align = "hv", 
          nrow = 3)

