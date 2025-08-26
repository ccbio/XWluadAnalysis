#BiocManager::install("sigminer", dependencies = TRUE)
library(sigminer)
library(NMF)
library(reshape2)
sbs96_sig_v3.2 <- as.data.frame(t(COSMIC_v3.2$signature$GRCh38$SBS96))
mut_transfer=function(mutl){
  x1=sub("^(.)(.)(.)(.)", "\\1[\\2>\\4]\\3", mutl)
  
  return(x1)
}
colnames(sbs96_sig_v3.2)=mut_transfer(colnames(sbs96_sig_v3.2))

sbs96_sig_v3.2_hg19 <- as.data.frame(t(COSMIC_v3.2$signature$GRCh37$SBS96))
colnames(sbs96_sig_v3.2_hg19)=mut_transfer(colnames(sbs96_sig_v3.2_hg19))

sig_mut=function(mut){
  top_rows_list <- list()
  for (col in 1:ncol(mut)) {
    current_col <- mut[, col]
    top_indices <- order(current_col, decreasing = TRUE)[1:16]
    top_row_names <- rownames(mut)[top_indices]
    top_rows_list[[colnames(mut)[col]]] <- top_row_names
  }
  return(top_rows_list)
}

#XW------------------------------------------------
mats <- mt_tally <- sig_tally(
  all_mut,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
  useSyn = TRUE,
  mode = "ALL"
)

str(mats, max.level = 1)
#mt_tally

#SBS-----------------------------------------------
#show_catalogue(mt_tally$SBS_96 %>% t(), mode = "SBS", style = "cosmic")
#est <- sig_estimate(mt_tally$SBS_96, range = 2:8, nrun = 2, verbose = TRUE)
#show_sig_number_survey2(est$survey)
#ggsave('fig1 clinical data/Fig1/signature/XW/NMF_select.pdf')
sigs <- sig_extract(mt_tally$SBS_96, n_sig = 4, nrun = 2)
p <- show_sig_profile(sigs, mode = "SBS", style = "cosmic")
p

signatures <- sigs$Signature.norm
weights <- sigs$Exposure.norm

# 计算每个特征的突变总数
total_mutations <- colSums(weights)

# 输出每个特征的突变总数
print(total_mutations)


#get_sig_similarity(sigs, sig_db = "legacy")
sim <- get_sig_similarity(sigs, sig_db = "SBS")
add_labels(p, x = 0.72, y = 0.25, y_end = 0.8, labels = sim, n_label = 4,font_size = 4)
ggsave('fig1 clinical data/Fig1/signature/XW/pattern_similarity.pdf')
sig_exposure_SBS_XW=get_sig_exposure(sigs)

show_sig_exposure(sigs, rm_space = TRUE, style = "cosmic")
ggsave('fig1 clinical data/Fig1/signature/XW/show_sig_exposure.pdf')
sigs_XW=sigs
df=sig_exposure_SBS_XW
df$sinature=NA
df$aetiology=NA
for (i in 1:nrow(df)){
  min_index <- which.max(df[i, 2:ncol(df)])
min_column_name <- names(df)[min_index+1]
df$aetiology[i]=sim$best_match[[min_column_name]]$aetiology
df$sinature[i]=strsplit(sim$best_match[[min_column_name]]$best_match,' ')[[1]][3]
}
sig_exposure_SBS_XW=df

sigs_Exposure.norm_XW=as.data.frame(sigs$Exposure.norm)
row.names(sigs_Exposure.norm_XW)=c(strsplit(sim$best_match[["Sig1"]]$best_match,' ')[[1]][3],
                                     strsplit(sim$best_match[["Sig2"]]$best_match,' ')[[1]][3],
                                     strsplit(sim$best_match[["Sig3"]]$best_match,' ')[[1]][3],
                                     strsplit(sim$best_match[["Sig4"]]$best_match,' ')[[1]][3])
sigs_Exposure.norm_XW_t=as.data.frame(t(sigs_Exposure.norm_XW))
sigs_Exposure.norm_XW_t$Tumor_Sample_Barcode=row.names(sigs_Exposure.norm_XW_t)
sigs_Exposure.norm_XW_t=melt(sigs_Exposure.norm_XW_t,id.vars = "Tumor_Sample_Barcode")

fig1_data=as.data.frame(mats$APOBEC_scores)
fig1_data_1=fig1_data[,c("Tumor_Sample_Barcode",'n_mutations')]
fig1_data_XW=merge(sigs_Exposure.norm_XW_t,fig1_data_1)

#sample_summary <- mafSummary(all_mut)
#sample_summary=sample_summary$variant.type.summary
#fig1_data_XW=merge(fig1_data_XW,sample_summary)
#####deconstructSigs----------------------------------
library(deconstructSigs)
mut_XW=read.csv('fig1 clinical data/Fig1/signature/hg38_all_mut.mut',sep='\t')
Signature.norm_mut=as.data.frame(sigs[["Signature.norm"]])
colnames(Signature.norm_mut)=c(strsplit(sim$best_match[["Sig1"]]$best_match,' ')[[1]][3],
                                   strsplit(sim$best_match[["Sig2"]]$best_match,' ')[[1]][3],
                                   strsplit(sim$best_match[["Sig3"]]$best_match,' ')[[1]][3],
                                   strsplit(sim$best_match[["Sig4"]]$best_match,' ')[[1]][3])

sig_mut_XW=sig_mut(Signature.norm_mut)
SBS4_XW=mut_XW[mut_XW$mut %in% sig_mut_XW$SBS4,]

whichSignatures_XW=data.frame()
sampl=row.names(mats$SBS_96)
sampl <- sampl[!sampl %in% c("1193067")]
for (i in sampl){
  x=whichSignatures(tumor.ref = as.data.frame(mats$SBS_96),
                signatures.ref = sbs96_sig_v3.2,  # 使用COSMIC数据库的特征
                sample.id = i,
                contexts.needed = TRUE,  # 是否需要上下文信息
                tri.counts.method = 'default')

  whichSignatures_XW=rbind(whichSignatures_XW,as.data.frame(x$weights))
}
whichSignatures_XW_clean <- whichSignatures_XW[, colMeans(whichSignatures_XW) > 0.01]


#EAS------------------------------------------------
test1=mixco_nosmokerF.maf[startsWith(mixco_nosmokerF.maf$Tumor_Sample_Barcode, "A"),]
test2=mixco_nosmokerF.maf[startsWith(mixco_nosmokerF.maf$Tumor_Sample_Barcode, "B"),]
TEST=rbind(test1,test2)

table(TEST$Chromosome)
write.table(TEST,'fig1 clinical data/Fig1/signature/EAS_signature.maf',sep='\t',quote = F,row.names = F)
EAS_nosmokerF_mut=read.maf(TEST)
mats <- mt_tally <- sig_tally(
  EAS_nosmokerF_mut,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  useSyn = TRUE,
  mode = "ALL"
)

str(mats, max.level = 1)
#mt_tally

#SBS-----------------------------------------------
#show_catalogue(mt_tally$SBS_96 %>% t(), mode = "SBS", style = "cosmic")
#est <- sig_estimate(mt_tally$SBS_96, range = 2:8, nrun = 2, verbose = TRUE)
#show_sig_number_survey2(est$survey)
#ggsave('fig1 clinical data/Fig1/signature/EAS/NMF_select.pdf')
sigs <- sig_extract(mt_tally$SBS_96, n_sig = 4, nrun = 2)
p <- show_sig_profile(sigs, mode = "SBS", style = "cosmic")
p

#get_sig_similarity(sigs, sig_db = "legacy")
sim <- get_sig_similarity(sigs, sig_db = "SBS")
add_labels(p, x = 0.72, y = 0.25, y_end = 0.8, labels = sim, n_label = 4,font_size = 4)
ggsave('fig1 clinical data/Fig1/signature/EAS/pattern_similarity.pdf')
sig_exposure_SBS_EAS=get_sig_exposure(sigs)

show_sig_exposure(sigs, rm_space = TRUE, style = "cosmic")
ggsave('fig1 clinical data/Fig1/signature/EAS/show_sig_exposure.pdf')
sigs_EAS=sigs
df=sig_exposure_SBS_EAS
df$sinature=NA
df$aetiology=NA
for (i in 1:nrow(df)){
  min_index <- which.max(df[i, 2:ncol(df)])
  min_column_name <- names(df)[min_index+1]
  df$aetiology[i]=sim$best_match[[min_column_name]]$aetiology
  df$sinature[i]=strsplit(sim$best_match[[min_column_name]]$best_match,' ')[[1]][3]
}
sig_exposure_SBS_EAS=df

sigs_Exposure.norm_EAS=as.data.frame(sigs$Exposure.norm)
row.names(sigs_Exposure.norm_EAS)=c(strsplit(sim$best_match[["Sig1"]]$best_match,' ')[[1]][3],
                                     strsplit(sim$best_match[["Sig2"]]$best_match,' ')[[1]][3],
                                     strsplit(sim$best_match[["Sig3"]]$best_match,' ')[[1]][3],
                                     strsplit(sim$best_match[["Sig4"]]$best_match,' ')[[1]][3])

sigs_Exposure.norm_EAS_t=as.data.frame(t(sigs_Exposure.norm_EAS))
sigs_Exposure.norm_EAS_t$Tumor_Sample_Barcode=row.names(sigs_Exposure.norm_EAS_t)
sigs_Exposure.norm_EAS_t=melt(sigs_Exposure.norm_EAS_t,id.vars = "Tumor_Sample_Barcode")

fig1_data=as.data.frame(mats$APOBEC_scores)
fig1_data_1=fig1_data[,c("Tumor_Sample_Barcode",'n_mutations')]
fig1_data_EAS=merge(sigs_Exposure.norm_EAS_t,fig1_data_1)

mut_EAS=read.csv('fig1 clinical data/Fig1/signature/EAS_signature.mut',sep='\t')
Signature.norm_mut=as.data.frame(sigs[["Signature.norm"]])
colnames(Signature.norm_mut)=c(strsplit(sim$best_match[["Sig1"]]$best_match,' ')[[1]][3],
                               strsplit(sim$best_match[["Sig2"]]$best_match,' ')[[1]][3],
                               strsplit(sim$best_match[["Sig3"]]$best_match,' ')[[1]][3],
                               strsplit(sim$best_match[["Sig4"]]$best_match,' ')[[1]][3])

sig_mut_EAS=sig_mut(Signature.norm_mut)
SBS4_EAS=mut_EAS[mut_EAS$mut %in% sig_mut_EAS$SBS4,]

whichSignatures_EAS=data.frame()
sampl=row.names(mats$SBS_96)
#sampl <- sampl[!sampl %in% c("1193067")]
for (i in sampl){
  x=whichSignatures(tumor.ref = as.data.frame(mats$SBS_96),
                    signatures.ref = sbs96_sig_v3.2_hg19,  # 使用COSMIC数据库的特征
                    sample.id = i,
                    contexts.needed = TRUE,  # 是否需要上下文信息
                    tri.counts.method = 'default')
  
  whichSignatures_EAS=rbind(whichSignatures_EAS,as.data.frame(x$weights))
}
whichSignatures_EAS_clean <- whichSignatures_EAS[, colMeans(whichSignatures_EAS) > 0.01]

#TCGA------------------------------------------------
TEST=mixco_nosmokerF.maf[startsWith(mixco_nosmokerF.maf$Tumor_Sample_Barcode, "TCGA"),]
write.table(TEST,'fig1 clinical data/Fig1/signature/TCGA_signature.maf',sep='\t',quote = F,row.names = F)
TCGA_nosmokerF_mut=read.maf(TEST)
mats <- mt_tally <- sig_tally(
  TCGA_nosmokerF_mut,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
  useSyn = TRUE,
  mode = "ALL"
)

str(mats, max.level = 1)
#mt_tally

#SBS-----------------------------------------------
show_catalogue(mt_tally$SBS_96 %>% t(), mode = "SBS", style = "cosmic")
#est <- sig_estimate(mt_tally$SBS_96, range = 2:8, nrun = 2, verbose = TRUE)
#show_sig_number_survey2(est$survey)
#ggsave('fig1 clinical data/Fig1/signature/TCGA/NMF_select.pdf')
sigs <- sig_extract(mt_tally$SBS_96, n_sig = 4, nrun = 2)
p <- show_sig_profile(sigs, mode = "SBS", style = "cosmic")
p

#get_sig_similarity(sigs, sig_db = "legacy")
sim <- get_sig_similarity(sigs, sig_db = "SBS")
add_labels(p, x = 0.72, y = 0.25, y_end = 0.8, labels = sim, n_label = 4,font_size = 4)
ggsave('fig1 clinical data/Fig1/signature/TCGA/pattern_similarity.pdf')
sig_exposure_SBS_TCGA=get_sig_exposure(sigs)

show_sig_exposure(sigs, rm_space = TRUE, style = "cosmic")
ggsave('fig1 clinical data/Fig1/signature/TCGA/show_sig_exposure.pdf')
sigs_TCGA=sigs
df=sig_exposure_SBS_TCGA
df$sinature=NA
df$aetiology=NA
for (i in 1:nrow(df)){
  min_index <- which.max(df[i, 2:ncol(df)])
  min_column_name <- names(df)[min_index+1]
  df$aetiology[i]=sim$best_match[[min_column_name]]$aetiology
  df$sinature[i]=strsplit(sim$best_match[[min_column_name]]$best_match,' ')[[1]][3]
}
sig_exposure_SBS_TCGA=df

sigs_Exposure.norm_TCGA=as.data.frame(sigs$Exposure.norm)
row.names(sigs_Exposure.norm_TCGA)=c(strsplit(sim$best_match[["Sig1"]]$best_match,' ')[[1]][3],
  strsplit(sim$best_match[["Sig2"]]$best_match,' ')[[1]][3],
  strsplit(sim$best_match[["Sig3"]]$best_match,' ')[[1]][3],
  strsplit(sim$best_match[["Sig4"]]$best_match,' ')[[1]][3])

sigs_Exposure.norm_TCGA_t=as.data.frame(t(sigs_Exposure.norm_TCGA))
sigs_Exposure.norm_TCGA_t$Tumor_Sample_Barcode=row.names(sigs_Exposure.norm_TCGA_t)
sigs_Exposure.norm_TCGA_t=melt(sigs_Exposure.norm_TCGA_t,id.vars = "Tumor_Sample_Barcode")

fig1_data=as.data.frame(mats$APOBEC_scores)
fig1_data_1=fig1_data[,c("Tumor_Sample_Barcode",'n_mutations')]
fig1_data_TCGA=merge(sigs_Exposure.norm_TCGA_t,fig1_data_1)

mut_TCGA=read.csv('fig1 clinical data/Fig1/signature/TCGA_signature.mut',sep='\t')
Signature.norm_mut=as.data.frame(sigs[["Signature.norm"]])
colnames(Signature.norm_mut)=c(strsplit(sim$best_match[["Sig1"]]$best_match,' ')[[1]][3],
                               strsplit(sim$best_match[["Sig2"]]$best_match,' ')[[1]][3],
                               strsplit(sim$best_match[["Sig3"]]$best_match,' ')[[1]][3],
                               strsplit(sim$best_match[["Sig4"]]$best_match,' ')[[1]][3])

sig_mut_TCGA=sig_mut(Signature.norm_mut)
SBS4_TCGA=mut_TCGA[mut_TCGA$mut %in% sig_mut_TCGA$SBS4,]

whichSignatures_TCGA=data.frame()
sampl=row.names(mats$SBS_96)
#sampl <- sampl[!sampl %in% c("1193067")]
for (i in sampl){
  x=whichSignatures(tumor.ref = as.data.frame(mats$SBS_96),
                    signatures.ref = sbs96_sig_v3.2,  # 使用COSMIC数据库的特征
                    sample.id = i,
                    contexts.needed = TRUE,  # 是否需要上下文信息
                    tri.counts.method = 'default')
  
  whichSignatures_TCGA=rbind(whichSignatures_TCGA,as.data.frame(x$weights))
}
whichSignatures_TCGA_clean <- whichSignatures_TCGA[, colMeans(whichSignatures_TCGA) >0.01]

#compare----------------------------------
##SBS cohort ========================================
library("pcutils")
signature_group=read.csv('fig1 clinical data/Fig1/signature/signature_group.txt',sep='\t')
row.names(signature_group)=signature_group$Signature.subgroup
sig_group=function(sig_group1){
  input_string <- signature_group[sig_group1,2]
  #print(input_string)
  # 按逗号拆分字符串
  split_result <- strsplit(input_string, ", ")[[1]]
  # 在每个元素前加上 "SBS"
  result_vector <- paste0("SBS", split_result)
  # 查看结果
  return(result_vector)
}
signature_group_SBS=list()
for (i in 1:nrow(signature_group)){
  signature_group_SBS[[signature_group[i,1]]]=sig_group(i)
}

whichsignature_gsva=function(whichsignature_matrix){
  for (i in names(signature_group_SBS)){
    #print(i)
    y=intersect(signature_group_SBS[[i]],colnames(whichsignature_matrix))
    if (length(y)==1){
      whichsignature_matrix[,i]=whichsignature_matrix[,y]
    }else{
      whichsignature_matrix[,i]=rowSums(whichsignature_matrix[,y])
    }
    
  }
  return(whichsignature_matrix)
}
whichsignature_gsva_XW=whichsignature_gsva(whichSignatures_XW)
whichsignature_gsva_EAS=whichsignature_gsva(whichSignatures_EAS)
whichsignature_gsva_TCGA=whichsignature_gsva(whichSignatures_TCGA)
whichsignature_gsva_XW$cohort='XW'
whichsignature_gsva_EAS$cohort='EAS'
whichsignature_gsva_TCGA$cohort='TCGA'

signature_gsva=rbind(whichsignature_gsva_XW,whichsignature_gsva_EAS,whichsignature_gsva_TCGA)
max_col_names <- apply(signature_gsva[,names(signature_group_SBS)], 
                       1, function(x) names(signature_gsva[,names(signature_group_SBS)])[which.max(x)])
signature_gsva$max_gsva=max_col_names
signature_gsva$signature=NA
ss=intersect(row.names(immune_matrix),row.names(signature_gsva))
signature_gsva[ss,'signature']=immune_matrix[match(ss,row.names(signature_gsva)),"signature" ]


plot_subtype_sig=function(dd,dd1){
  meta1=signature_gsva
  #meta1=signature_gsva[signature_gsva$cluster %in% aa,]
  meta1=meta1[,c("cohort", dd)]
  meta1=na.omit(meta1)
  max_ = max(meta1[,dd]) #+ max(tmb_all$total_perMB)/3 
  #meta1$cluster=as.factor(gsub(paste0(bb,'_'),'',meta1$cluster))
  #print(max_)
  df1_g <- as.character(unique(meta1$cohort))
  # 使用内置函数combn(x, m, FUN, simplify)获取比较组
  # x 要取组合的对象，为一个vector
  # m 要取出的数量
  # FUN为对取出的每个组合之行的函数，默认为NULL，即不执行
  # simplify为是否简单化输出，默认为T，输出的是data frame；若为F，输出的是list
  df1_cmp <- combn(df1_g, 2, simplify = F) 
  meta1$cluster=factor(meta1$cohort,levels = c('XW','EAS','TCGA'))
  p <- ggviolin(meta1, "cohort", dd , 
                fill = "cohort", #小提琴内部颜色对应的数据列
                color = "cohort", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
                trim = T, 
                palette = c('#E5A094','#148486','#EEB033'),#c("#015493","#019092","#999999","#F4A99B"),  #小提琴自定义颜色
                legend = "right", #图例添加在图的右侧
                legend.title = " ",#图例的标题
                font.legend = c(12, "plain", "black"), #图例字体的大小/样式/颜色
                font.y = 12,  #y轴标题字体大小
                
                size=0.25,
                font.tickslab = c(12,"plain","black"), #x轴 y轴刻度大小/样式/颜色
                add = "boxplot",  #叠加箱线图
                add.params = list(   
                  fill = "white", #设置箱线图内部颜色
                  color = "black",  #设置箱线图边框
                  width = 0.1,   #箱线图的宽度
                  linetype = 0.5,
                  size=0.25)) +
    labs(x = NULL, #设置x轴标题
         y = dd1, #设置y轴标题
         title = dd1) + 
    #geom_signif(  #添加显著性标记
    #  comparisons = df1_cmp ,#指定比较对象
    #  map_signif_level = T ,  #指定显著性差异表示方式，布尔型变量，如果为TRUE，就用***形式来展示显著性差异，"***"=0.001, "**"=0.01, "*"=0.05
    #  textsize = 3 , #标记文字的大小  
    #  test = wilcox.test, #指定检验方法  含wilcox.test、 t.test
    #  step_increase = 0.05 , #多个组时，差异标注的距离 
    #  size = 1 , #标记线条的粗线
    #  tip_length = 0.02,   #标记短竖线的长度
    #  y_position = max_#标记在纵轴方向上的位置
    #)+
    geom_signif(
      comparisons = df1_cmp,
      test = "wilcox.test", # 你可以选择适合你数据的测试
      map_signif_level = function(p) {
        
        if (p<0.001){
          return('***')}
        else if (p<0.01){
          return('**')}
        else if (p<0.05){
          return('*')
        }
        
        else {
          return(NA) # 不显著时返回NA
        }
      },
      textsize = 5 , #标记文字的大小  
      step_increase = 0.05 , #多个组时，差异标注的距离 
      size = 0.25 , #标记线条的粗线
      tip_length = 0.02,   #标记短竖线的长度
      y_position = max_,#标记在纵轴方向上的位置
      vjust = 0.5
    )+
    theme(legend.position='none',
          axis.line = element_line(size=0.25),
          axis.title.y = element_blank(),
          panel.border = element_rect(color = "black", size = 0.25, fill = NA),
          axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75),
          plot.title=element_text(size=10,
                                  face="bold",
                                  color="black", #颜色
                                  hjust=0.5, #调整位置，正中间
                                  lineheight=0.5)
          #plot.background = element_rect(color = "blue", size = 2)
    )
  return(p)
}

pl=list()
for (i in names(signature_group_SBS)[-c(11,13)]){
  print(i)
  pl[[i]]=plot_subtype_sig(i,i)
}
P2=grid.arrange(grobs = pl,ncol = 6, nrow = 2)

pl1=list()
pl1[['ALL']]=ggscatter(data=signature_gsva[signature_gsva$signature %in% c(SBS[[c("Tobacco smoking")]],SBS[["APOBEC3 activity"]]),],x="Tobacco_signatures",
                       y="APOBEC_signatures",color = "signature",shape = "cohort",
                       size=3,
                       palette = color_SG,
                       ggtheme=theme_minimal()
)+labs(#x = NULL, #设置x轴标题
  #y = dd1, #设置y轴标题
  title = 'Three Cohorts') + 
  stat_cor(method="pearson",
           label.x=0.2,label.y=0.5)+
  geom_smooth(aes(x = Tobacco_signatures, y = APOBEC_signatures, group = 1), 
              method = "lm", se = FALSE, color = "blue")+
  theme(legend.position='none',
        axis.line = element_line(size=0.25),
        #axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 13),
        plot.title=element_text(size=13,
                                face="bold",
                                color="black", #颜色
                                hjust=0.5, #调整位置，正中间
                                lineheight=0.5)
        #plot.background = element_rect(color = "blue", size = 2)
  )
for (i in c('XW','EAS','TCGA')){
  pl1[[i]]=ggscatter(data=signature_gsva[signature_gsva$cohort==i&signature_gsva$signature %in% c(SBS[[c("Tobacco smoking")]],SBS[["APOBEC3 activity"]]),],x="Tobacco_signatures",
          y="APOBEC_signatures",color = "signature",shape = "cohort",palette = color_SG,
          size=3,
          ggtheme=theme_minimal()
)+labs(#x = NULL, #设置x轴标题
       #y = dd1, #设置y轴标题
       title = i) + 
  stat_cor(method="pearson",
           label.x=0.2,label.y=0.5)+
  geom_smooth(aes(x = Tobacco_signatures, y = APOBEC_signatures, group = 1), 
              method = "lm", se = FALSE, color = "blue")+
    theme(legend.position='none',
          axis.line = element_line(size=0.25),
          #axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", size = 0.25, fill = NA),
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(color = "black", size = 13),
          plot.title=element_text(size=13,
                                  face="bold",
                                  color="black", #颜色
                                  hjust=0.5, #调整位置，正中间
                                  lineheight=0.5)
          #plot.background = element_rect(color = "blue", size = 2)
    )
}

P3=grid.arrange(grobs = pl1,ncol = 1, nrow = 4)
##SBS4_GENE======================================================
aa=list(EAS=intersect(unique(SBS4_EAS$Hugo_Symbol),gene_XW),
XW=intersect(unique(SBS4_XW$Hugo_Symbol),gene_XW),
TCGA=intersect(unique(SBS4_TCGA$Hugo_Symbol),gene_XW))


##属于某个SBS的人的drivergene 突变数量=================================
driver_data=function(data){
  #samp=c()
  for (gg in gene_XW){
    samp=c(unlist(unique(all_mut@data[all_mut@data$Hugo_Symbol==gg,
                           "Tumor_Sample_Barcode"])),
         unlist(unique(TCGA_nosmokerF_mut@data[TCGA_nosmokerF_mut@data$Hugo_Symbol==gg,
                               "Tumor_Sample_Barcode"])),
         unlist(unique(EAS_nosmokerF_mut@data[EAS_nosmokerF_mut@data$Hugo_Symbol==gg,
                               "Tumor_Sample_Barcode"])))
  
  
  
  data[,gg]=0
  data[data$Tumor_Sample_Barcode %in% samp,gg]=1}
  return(data)
}
drivers_sig=drivers_sig=driver_data(immune_matrix)
drivers_sig$per=1
drivers_sig$signature=factor(drivers_sig$signature,levels = c("SBS4","SBS2","SBS13","SBS3","SBS6","SBS1","SBS5"))
drivers_sig$Cohort=factor(drivers_sig$Cohort,levels = c('XW','EAS','TCGA'))

plot_driver_sig=function(gg){
  x=paste0(gg," mut type") 
  df=drivers_sig[drivers_sig[,gg]==1,]
  df$Cohort=factor(df$Cohort,levels = c('XW','EAS','TCGA'))
  df_summary <-  df%>%
    group_by(Cohort) %>%
    summarise(N = n()) %>%
    as.data.frame
  df$cc=NA
  for (i in 1:nrow(df_summary)){
    df[df$Cohort==df_summary[i,1],'cc']=paste0(df_summary[i,1],'\n',
                                               '(n = ',df_summary[i,2],')')
  }
  
  df$cc=factor(df$cc,levels = unique(df$cc))
  p=ggplot(df, aes(x = cc, y = per, fill = signature)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = color_SG) +
    labs(x = "", 
        y = "Proportion of SBS signature", 
         fill = "Signature",
        title = x) +
    theme_bw()+
    theme(#legend.position='none',
          panel.grid= element_blank(),
          axis.line = element_line(size=0.25),
          axis.title = element_text(color = "black",face="bold",size = 13),
          #panel.border = element_rect(color = "black", size = 0.25, fill = NA),
          axis.text = element_text(color = "black",size = 12),
          plot.title=element_text(size=15,
                                  face="bold",
                                  color="black", #颜色
                                  hjust=0.5, #调整位置，正中间
                                  lineheight=0.5)
          #plot.background = element_rect(color = "blue", size = 2)
  )
  
return(p)
}
pl=list()
#for (gg in c("TP53","RYR2","KRAS",'EGFR')){
for (gg in gene_XW){
  print(gg)
  pl[[gg]]=plot_driver_sig(gg)
}
P4=grid.arrange(grobs = pl,ncol = 6, nrow = 8)

pl=list()
for (i in c("APOBEC_signatures","Tobacco_signatures")){
  print(i)
  pl[[i]]=plot_subtype_sig(i,i)
}

for (gg in c("TP53","RYR2","KRAS",'EGFR')){
  print(gg)
  pl[[gg]]=plot_driver_sig(gg)
}
P4=grid.arrange(grobs = pl,ncol = 4, nrow = 1)



##mutation and expose===============================
fig1_data_TCGA$cohort='TCGA'
fig1_data_EAS$cohort='EAS'
fig1_data_XW$cohort='XW'

fig1_data=rbind(fig1_data_XW,fig1_data_EAS,fig1_data_TCGA)
fig1_data$n_mutations=log10(fig1_data$n_mutations)
fig1_data$variable=factor(fig1_data$variable,level=c('SBS4',
                                                     'SBS2','SBS13',
                                                     'SBS3','SBS6',
                                                     'SBS1','SBS5'))
p1=ggscatter(fig1_data[fig1_data$value>0.0001,],x="value",y="n_mutations",color = "variable",shape = "cohort",
          palette = color_SG,
          ggtheme=theme_minimal(),
          size = 2
)+labs(x='Activity',y='Number of mutations (log10)')+
  stat_cor(method="pearson",
           label.x=0.5,label.y=4)+
  geom_smooth(aes(x = value, y = n_mutations, group = 1), 
              method = "lm", se = FALSE, color = "blue")+
  theme_bw()+
  theme(#legend.position='none',
        axis.line = element_line(size=0.25),
        axis.title = element_text(size=12,
                                    face="bold",
                                    color="black"),
        panel.grid.major = element_blank(),   # 删除主要网格线
        panel.grid.minor = element_blank(),   # 删除次要网格线
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        axis.text = element_text(size=12,color = "black"),
        plot.title=element_text(size=12,
                                face="bold",
                                color="black", #颜色
                                hjust=0.5, #调整位置，正中间
                                lineheight=0.5)
        #plot.background = element_rect(color = "blue", size = 2)
  )


  max_ = max(fig1_data$n_mutations) #+ max(tmb_all$total_perMB)/3 
  df1_g <- as.character(unique(fig1_data$variable))
  df1_cmp <- combn(df1_g, 2, simplify = F) 
fig1_data$cohort=factor(fig1_data$cohort,levels = c('XW','EAS','TCGA'))
pl=list()
SBS=list()
SBS[['Tobacco smoking']]=c("SBS4")
SBS[['APOBEC3 activity']]=c("SBS2","SBS13")
SBS[['DNA repair']]=c("SBS3","SBS6")
SBS[['Clock−like']]=c("SBS1","SBS5")
for (i in names(SBS)){
  print(i)
  pl[[i]]=ggscatter(fig1_data[fig1_data$value>0.0001&fig1_data$variable %in% SBS[[i]],],x="value",y="n_mutations",color = "variable",shape = "cohort",
                    palette = color_SG,
                    ggtheme=theme_minimal(),
                    size = 2
  )+labs(x='Activity',y='Number of mutations (log10)',title = i)+
    stat_cor(method="pearson",
             label.x=0.5,label.y=4)+
    geom_smooth(aes(x = value, y = n_mutations, group = 1), 
                method = "lm", se = FALSE, color = "blue")+
    theme_bw()+
    theme(#legend.position='none',
      axis.line = element_line(size=0.25),
      axis.title = element_text(size=12,
                                face="bold",
                                color="black"),
      panel.grid.major = element_blank(),   # 删除主要网格线
      panel.grid.minor = element_blank(),   # 删除次要网格线
      panel.border = element_rect(color = "black", size = 0.25, fill = NA),
      axis.text = element_text(size=12,color = "black"),
      plot.title=element_text(size=12,
                              face="bold",
                              color="black", #颜色
                              hjust=0.5, #调整位置，正中间
                              lineheight=0.5)
      #plot.background = element_rect(color = "blue", size = 2)
    )
  
}
P2=grid.arrange(grobs = pl,ncol = 2, nrow = 2)


p2=ggviolin(fig1_data, "variable", "n_mutations" , fill="variable",
           trim = T, 
           palette = color_SG,#c("#015493","#019092","#999999","#F4A99B"),  #小提琴自定义颜色
           legend = "right", #图例添加在图的右侧
           legend.title = " ",#图例的标题
           font.legend = c(8, "plain", "black"), #图例字体的大小/样式/颜色
           font.y = 8,  #y轴标题字体大小
           
           size=0.25,alpha = 0.8,
           font.tickslab = c(8,"plain","black"),
         add = "boxplot",  #叠加箱线图
         add.params = list(   
           fill = "white", #设置箱线图内部颜色
           color = "black",  #设置箱线图边框
           width = 0.25,   #箱线图的宽度
           linetype = 1,
           size=0.25))+
  labs(y = 'Number of mutations (log10)', #设置x轴标题
       x = 'Signatures') + 
  geom_signif(
    comparisons = df1_cmp,
    test = "wilcox.test", # 你可以选择适合你数据的测试
    map_signif_level = function(p) {
      
      if (p<0.001){
        return('***')}
      else if (p<0.01){
        return('**')}
      else if (p<0.05){
        return('*')
      }
      
      else {
        return(NA) # 不显著时返回NA
      }
    },
    textsize = 5 , #标记文字的大小  
    step_increase = 0.05 , #多个组时，差异标注的距离 
    size = 0.25 , #标记线条的粗线
    tip_length = 0.02,   #标记短竖线的长度
    y_position = max_,#标记在纵轴方向上的位置
    vjust = 0.5
  )+
  geom_point(aes(x = variable, y = n_mutations, color = cohort),fill='black',
    position = position_jitterdodge(jitter.width = 0.1), size = 2,alpha=0.8) +  # 添加分组散点图
  scale_color_manual(values = c("#D07459","#E2BF6B","#2E968A")) +
  theme_bw()+
  theme(#legend.position='none',
    axis.line = element_line(size=0.25),
    axis.title = element_text(size=12,
                              face="bold",
                              color="black"),
    panel.grid.major = element_blank(),   # 删除主要网格线
    panel.grid.minor = element_blank(),   # 删除次要网格线
    panel.border = element_rect(color = "black", size = 0.25, fill = NA),
    axis.text = element_text(size=12,color = "black"),
    plot.title=element_text(size=12,
                            face="bold",
                            color="black", #颜色
                            hjust=0.5, #调整位置，正中间
                            lineheight=0.5)
    #plot.background = element_rect(color = "blue", size = 2)
  )

#pie================================
sig_exposure_SBS_TCGA$cohort='TCGA'
sig_exposure_SBS_EAS$cohort='EAS'
sig_exposure_SBS_XW$cohort='XW'

sig_exposure_SBS=rbind(sig_exposure_SBS_TCGA,sig_exposure_SBS_EAS,sig_exposure_SBS_XW)

sig_exposure_SBS$SG=NA
sig_exposure_SBS$SG[grepl("APOBEC",sig_exposure_SBS$aetiology)]='APOBEC'
sig_exposure_SBS$SG[grepl("DNA",sig_exposure_SBS$aetiology)]='DNA repair'
sig_exposure_SBS$SG[grepl("smoking",sig_exposure_SBS$aetiology)]='Tobacco smoking'
sig_exposure_SBS$SG[grepl("clock-like",sig_exposure_SBS$aetiology)]='clock-like signature'

color_signature=c("#F8CDCF","#C2E99E","#FAC5A4","#FBEAA6","#ABEAF9",'#AAA9F7')
names(color_signature)=unique(sig_exposure_SBS$aetiology)

pie_sig=function(dd,aa,bb){

  df <- as.data.frame(table(dd$sinature)) %>% 
    magrittr::set_colnames( c("class", "freq")) %>% 
  mutate(per = freq/sum(freq)) %>% 
  mutate(label = paste0(round(per,4)*100,'%'))
#df
  df$class=factor(df$class,levels = bb)
  p1 = ggplot(df, aes(x = 2, y=freq, fill = class)) + 
  geom_bar(width = 1, stat = "identity")+
  xlim(1,2.5)+
  geom_text(aes(label = paste0(class,': ',label)),position =  position_stack(vjust = 0.5))+
  scale_fill_manual(values = color_SG) +
  theme_classic(base_size = 12)+
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5)) + 
  coord_polar(theta = "y", start=0)+
  labs(fill="class",  x=NULL,  y=NULL, 
       title=paste0(aa,'(n=',sum(df$freq),")"))+
  theme(legend.position="none")

  return(p1)
  }
#"Tobacco smoking"（SBS4）
#"Activity of APOBEC family of cytidine deaminases"（SBS13、SBS2）
#DNA repair （SBS3、SBS6）
#clock-like signature(SBS1、SBS5)
p1=pie_sig(sig_exposure_SBS_XW,'XW',c('SBS4','SBS2','SBS5','SBS6'))
p2=pie_sig(sig_exposure_SBS_EAS,'EAS',c('SBS4','SBS2','SBS1','SBS5'))
p3=pie_sig(sig_exposure_SBS_TCGA,'TCGA',c('SBS4','SBS13','SBS3','SBS6'))


p1|p2|p3

df=sigs_Exposure.norm_TCGA
df$signature=row.names(df)
data1 <- melt(df, id = 'signature')
data1$cohort='TCGA'

df=sigs_Exposure.norm_XW
df$signature=row.names(df)
data2 <- melt(df, id = 'signature')
data2$cohort='XW'

df=sigs_Exposure.norm_EAS
df$signature=row.names(df)
data3 <- melt(df, id = 'signature')
data3$cohort='EAS'

sigs_Exposure.norm=rbind(data2,data3,data1)
#"Tobacco smoking"（SBS4）
#"Activity of APOBEC family of cytidine deaminases"（SBS13、SBS2）
#DNA repair （SBS3、SBS6）
#clock-like signature(SBS1、SBS5) 

color_SG=c("#BD554E","#FFF3D1","#88B37E","#AEDFE4","#AEDFE4","#88B37E","#FFF3D1")

names(color_SG)=unique(sigs_Exposure.norm$signature)

sigs_Exposure.norm$signature=factor(sigs_Exposure.norm$signature,
                                    levels = c('SBS4','SBS2','SBS13','SBS1','SBS5','SBS3','SBS6'))
sigs_Exposure.norm$cohort=factor(sigs_Exposure.norm$cohort,
                                    levels = c('XW','EAS','TCGA'))

p5=ggplot(sigs_Exposure.norm,aes(x = variable, y = value))+  
  geom_col(aes(fill = signature),width = 1)+  
  scale_fill_manual(values = color_SG)+
  theme_minimal() +
  theme_prism(base_fontface = "plain",              
              base_family = "serif",               
              base_size = 16,                
              base_line_size = 0.8,               
              axis_text_angle = 0) +  # 取消背景板网格线，并添加图片框线  
  theme(panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),  
        axis.ticks.x = element_blank(),
        #axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),        
        axis.title.y = element_text(face="bold",size = 14, 
                                    color = "black",family = "sans"))#修改Y轴标题字体、大小、加p


P1=(p1|p2|p3)/p5


#ggviolin-------------------------------------------
cohort4_mixcopdata_Nofemale$signature=NA

row.names(cohort4_mixcopdata_Nofemale)=cohort4_mixcopdata_Nofemale$Tumor_Sample_Barcode
X=intersect(sig_exposure_SBS_XW$sample,cohort4_mixcopdata_Nofemale$Tumor_Sample_Barcode)
cohort4_mixcopdata_Nofemale[X,'signature']=sig_exposure_SBS_XW[match(X,sig_exposure_SBS_XW$sample),'sinature']
X=intersect(sig_exposure_SBS_EAS$sample,cohort4_mixcopdata_Nofemale$Tumor_Sample_Barcode)
cohort4_mixcopdata_Nofemale[X,'signature']=sig_exposure_SBS_EAS[match(X,sig_exposure_SBS_EAS$sample),'sinature']
X=intersect(sig_exposure_SBS_TCGA$sample,cohort4_mixcopdata_Nofemale$Tumor_Sample_Barcode)
cohort4_mixcopdata_Nofemale[X,'signature']=sig_exposure_SBS_TCGA[match(X,sig_exposure_SBS_TCGA$sample),'sinature']
cohort4_mixcopdata_Nofemale_signature=cohort4_mixcopdata_Nofemale[,1:10]
cohort4_mixcopdata_Nofemale_signature=na.omit(cohort4_mixcopdata_Nofemale_signature)

cohort4_mixcopdata_Nofemale$TMB=NA
X=intersect(tmb_all$Tumor_Sample_Barcode,cohort4_mixcopdata_Nofemale$Tumor_Sample_Barcode)
cohort4_mixcopdata_Nofemale[X,'TMB']=tmb_all[match(X,tmb_all$Tumor_Sample_Barcode),'total_perMB']

cohort4_mixcopdata_Nofemale$GII=NA
X=intersect(GII$Tumor_Sample_Barcode,cohort4_mixcopdata_Nofemale$Tumor_Sample_Barcode)
cohort4_mixcopdata_Nofemale[X,'GII']=GII[match(X,GII$Tumor_Sample_Barcode),'GII']

cohort4_mixcopdata_Nofemale_signature=cohort4_mixcopdata_Nofemale[,c(1:10,13,14)]
cohort4_mixcopdata_Nofemale_signature=na.omit(cohort4_mixcopdata_Nofemale_signature)

sig_cil=function(df,aa,B){
  #df=cohort4_mixcopdata_Nofemale_signature
df1_g <- as.character(unique(df$Cohort))
max_ = max(df$Age) #+ max(tmb_all$total_perMB)/3 
# 使用内置函数combn(x, m, FUN, simplify)获取比较组
# x 要取组合的对象，为一个vector
# m 要取出的数量
# FUN为对取出的每个组合之行的函数，默认为NULL，即不执行
# simplify为是否简单化输出，默认为T，输出的是data frame；若为F，输出的是list
df1_cmp <- combn(df1_g, 2, simplify = F) 
ann_colors_mixco$Cohort=c('#D07358',"#E2BF6A","#2B9789")
names(ann_colors_mixco$Cohort)=c('XW','EAS','TCGA')
p <- ggviolin(df, "Cohort", B , 
              fill = "Cohort", #小提琴内部颜色对应的数据列
              color = "Cohort", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
              trim = T, 
              palette = ann_colors_mixco$Cohort,  #小提琴自定义颜色
              legend = "right", #图例添加在图的右侧
              legend.title = " ",#图例的标题
              font.legend = c(15, "plain", "black"), #图例字体的大小/样式/颜色
              font.y = 15,  #y轴标题字体大小
              font.tickslab = c(15,"plain","black"), #x轴 y轴刻度大小/样式/颜色
              add = "boxplot",  #叠加箱线图
              add.params = list(   
                fill = "white", #设置箱线图内部颜色
                color = "black",  #设置箱线图边框
                width = 0.05,   #箱线图的宽度
                linetype = 1)) +
  labs(x = NULL, #设置x轴标题
       y = "Age (years)", #设置y轴标题
       title = aa) + 
  geom_signif(  #添加显著性标记
    comparisons = df1_cmp ,#指定比较对象
    map_signif_level = F ,  #指定显著性差异表示方式，布尔型变量，如果为TRUE，就用***形式来展示显著性差异，"***"=0.001, "**"=0.01, "*"=0.05
    textsize = 5 , #标记文字的大小  
    test = t.test, #指定检验方法  含wilcox.test、 t.test
    step_increase = 0.1 , #多个组时，差异标注的距离 
    size = 1 , #标记线条的粗线
    tip_length = 0.02,   #标记短竖线的长度
    y_position = max_#标记在纵轴方向上的位置
  )+
  theme(
    axis.line = element_line(size=1),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    #plot.background = element_rect(color = "blue", size = 2)
  )+theme(legend.position="none")
return(p)

}

#"Tobacco smoking"（SBS4）
#"Activity of APOBEC family of cytidine deaminases"（SBS13、SBS2）
#DNA repair （SBS3、SBS6）
#clock-like signature(SBS1、SBS5) 
p1=sig_cil(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c('SBS4'),],"Tobacco smoking","Age")
p2=sig_cil(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS13","SBS2"),],"APOBEC3 activity","Age")
p3=sig_cil(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS3","SBS6"),],"DNA repair","Age")
p4=sig_cil(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS1","SBS5"),],"clock-like signature","Age")
P2=(p1|p2)|(p3|p4)

sig_cil1=function(df,aa,B){
  #df=cohort4_mixcopdata_Nofemale_signature
  df1_g <- as.character(unique(df$Cohort))
  max_ = max(df$TMB) #+ max(tmb_all$total_perMB)/3 
  # 使用内置函数combn(x, m, FUN, simplify)获取比较组
  # x 要取组合的对象，为一个vector
  # m 要取出的数量
  # FUN为对取出的每个组合之行的函数，默认为NULL，即不执行
  # simplify为是否简单化输出，默认为T，输出的是data frame；若为F，输出的是list
  df1_cmp <- combn(df1_g, 2, simplify = F) 
  ann_colors_mixco$Cohort=c('#D07358',"#E2BF6A","#2B9789")
  names(ann_colors_mixco$Cohort)=c('XW','EAS','TCGA')
  p <- ggviolin(df, "Cohort", B , 
                fill = "Cohort", #小提琴内部颜色对应的数据列
                color = "Cohort", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
                trim = T, 
                palette = ann_colors_mixco$Cohort,  #小提琴自定义颜色
                legend = "right", #图例添加在图的右侧
                legend.title = " ",#图例的标题
                font.legend = c(15, "plain", "black"), #图例字体的大小/样式/颜色
                font.y = 15,  #y轴标题字体大小
                font.tickslab = c(15,"plain","black"), #x轴 y轴刻度大小/样式/颜色
                add = "boxplot",  #叠加箱线图
                add.params = list(   
                  fill = "white", #设置箱线图内部颜色
                  color = "black",  #设置箱线图边框
                  width = 0.05,   #箱线图的宽度
                  linetype = 1)) +
    labs(x = NULL, #设置x轴标题
         y = "TMB per MB", #设置y轴标题
         title = aa) + 
    geom_signif(  #添加显著性标记
      comparisons = df1_cmp ,#指定比较对象
      map_signif_level = F ,  #指定显著性差异表示方式，布尔型变量，如果为TRUE，就用***形式来展示显著性差异，"***"=0.001, "**"=0.01, "*"=0.05
      textsize = 5 , #标记文字的大小  
      test = t.test, #指定检验方法  含wilcox.test、 t.test
      step_increase = 0.1 , #多个组时，差异标注的距离 
      size = 1 , #标记线条的粗线
      tip_length = 0.02,   #标记短竖线的长度
      y_position = max_#标记在纵轴方向上的位置
    )+
    theme(
      axis.line = element_line(size=1),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      #plot.background = element_rect(color = "blue", size = 2)
    )+
    theme(legend.position="none")
  return(p)
  
} 
p1=sig_cil1(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c('SBS4'),],"Tobacco smoking","TMB")
p2=sig_cil1(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS13","SBS2"),],"APOBEC3 activity","TMB")
p3=sig_cil1(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS3","SBS6"),],"DNA repair","TMB")
p4=sig_cil1(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS1","SBS5"),],"clock-like signature","TMB")
P3=(p1|p2)|(p3|p4)
P3

sig_cil2=function(df,aa,B){
  #df=cohort4_mixcopdata_Nofemale_signature
  df1_g <- as.character(unique(df$Cohort))
  max_ = max(df$GII) #+ max(tmb_all$total_perMB)/3 
  # 使用内置函数combn(x, m, FUN, simplify)获取比较组
  # x 要取组合的对象，为一个vector
  # m 要取出的数量
  # FUN为对取出的每个组合之行的函数，默认为NULL，即不执行
  # simplify为是否简单化输出，默认为T，输出的是data frame；若为F，输出的是list
  df1_cmp <- combn(df1_g, 2, simplify = F) 
  ann_colors_mixco$Cohort=c('#D07358',"#E2BF6A","#2B9789")
  names(ann_colors_mixco$Cohort)=c('XW','EAS','TCGA')
  p <- ggviolin(df, "Cohort", B , 
                fill = "Cohort", #小提琴内部颜色对应的数据列
                color = "Cohort", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
                trim = T, 
                palette = ann_colors_mixco$Cohort,  #小提琴自定义颜色
                legend = "right", #图例添加在图的右侧
                legend.title = " ",#图例的标题
                font.legend = c(15, "plain", "black"), #图例字体的大小/样式/颜色
                font.y = 15,  #y轴标题字体大小
                font.tickslab = c(15,"plain","black"), #x轴 y轴刻度大小/样式/颜色
                add = "boxplot",  #叠加箱线图
                add.params = list(   
                  fill = "white", #设置箱线图内部颜色
                  color = "black",  #设置箱线图边框
                  width = 0.05,   #箱线图的宽度
                  linetype = 1)) +
    labs(x = NULL, #设置x轴标题
         y = "GII", #设置y轴标题
         title = aa) + 
    geom_signif(  #添加显著性标记
      comparisons = df1_cmp ,#指定比较对象
      map_signif_level = F ,  #指定显著性差异表示方式，布尔型变量，如果为TRUE，就用***形式来展示显著性差异，"***"=0.001, "**"=0.01, "*"=0.05
      textsize = 5 , #标记文字的大小  
      test = t.test, #指定检验方法  含wilcox.test、 t.test
      step_increase = 0.1 , #多个组时，差异标注的距离 
      size = 1 , #标记线条的粗线
      tip_length = 0.02,   #标记短竖线的长度
      y_position = max_#标记在纵轴方向上的位置
    )+
    theme(
      axis.line = element_line(size=1),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      #plot.background = element_rect(color = "blue", size = 2)
    )+
    theme(legend.position="none")
  return(p)
  
} 
p1=sig_cil2(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c('SBS4'),],"Tobacco smoking","GII")
p2=sig_cil2(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS13","SBS2"),],"APOBEC3 activity","GII")
p3=sig_cil2(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS3","SBS6"),],"DNA repair","GII")
p4=sig_cil2(cohort4_mixcopdata_Nofemale_signature[cohort4_mixcopdata_Nofemale_signature$signature %in% c("SBS1","SBS5"),],"clock-like signature","GII")
P4=(p1|p2)|(p3|p4)
P4

P1|(P2/P3/P4)
#mutation type========================================================



EAS_mut_Nofemale[,1:44]=lapply(EAS_mut_Nofemale[,1:44], as.numeric)
EAS_mut_Nofemale$Mut_count_in_driver=rowSums(EAS_mut_Nofemale[,1:44])
TCGA_mut_Nofemale[,1:44]=lapply(TCGA_mut_Nofemale[,1:44], as.numeric)
TCGA_mut_Nofemale$Mut_count_in_driver=rowSums(TCGA_mut_Nofemale[,1:44])
XW_mut_Nofemale[,1:44]=lapply(XW_mut_Nofemale[,1:44], as.numeric)
XW_mut_Nofemale$Mut_count_in_driver=rowSums(XW_mut_Nofemale[,1:44])
mut_Nofemale=rbind(EAS_mut_Nofemale,TCGA_mut_Nofemale,XW_mut_Nofemale)

immune_matrix$Mut_count_in_driver=NA
ss=intersect(immune_matrix$Tumor_Sample_Barcode,mut_Nofemale$Tumor_Sample_Barcode)
#pdataxwfa_new2[ss,"gmm_cluster"]=x[ss,"drmod$classification" ]
immune_matrix[ss,'Mut_count_in_driver']=mut_Nofemale[match(ss,mut_Nofemale$Tumor_Sample_Barcode),'Mut_count_in_driver']

SBS4_immune_matrix=immune_matrix[immune_matrix$signature %in% c('SBS4'),]

feature='Mut_count_in_driver'
x1=immune_matrix[,c('Age',feature,'signature')]
x1=na.omit(x1)
x1=x1[x1$signature=='SBS4'&x1$Age<=65,]
ggscatter(x1,x='Age',y=feature,
          color = 'signature',
            palette = color_SG,
            ggtheme=theme_minimal(),
            size = 2
  )+labs(x='Age',y='Number of mutations in driver genes',title = i)+
  stat_cor(method="pearson",
           label.x=50,label.y=15)+
  geom_smooth(aes(x=Age,y=x1[,feature], group = 1), 
              method = "lm", se = FALSE, color = "blue")+
  theme_bw()+
  theme(#legend.position='none',
    axis.line = element_line(size=0.25),
    axis.title = element_text(size=12,
                              face="bold",
                              color="black"),
    panel.grid.major = element_blank(),   # 删除主要网格线
    panel.grid.minor = element_blank(),   # 删除次要网格线
    panel.border = element_rect(color = "black", size = 0.25, fill = NA),
    axis.text = element_text(size=12,color = "black"),
    plot.title=element_text(size=12,
                            face="bold",
                            color="black", #颜色
                            hjust=0.5, #调整位置，正中间
                            lineheight=0.5)
    #plot.background = element_rect(color = "blue", size = 2)
  )



anno_maf_EAS=read.maf(read.csv('anno_maf/EAS_BARCODE.maf',sep='\t'))
anno_maf_TCGA=read.maf(read.csv('anno_maf/TCGA_nosmokerF_mut_BARCODE.maf',sep='\t'))

XW_HGVSp=as.data.frame(t(read.csv('anno_maf/XW_HGVSp.csv',check.names = F)))
colnames(XW_HGVSp)=XW_HGVSp[1,]
XW_HGVSp=XW_HGVSp[-1,]
XW_HGVSp=XW_HGVSp[,intersect(gene_XW,colnames(XW_HGVSp))]
EAS_HGVSp=as.data.frame(t(read.csv('anno_maf/EAS_HGVSp.csv',check.names = F)))
colnames(EAS_HGVSp)=EAS_HGVSp[1,]
EAS_HGVSp=EAS_HGVSp[-1,]
EAS_HGVSp=EAS_HGVSp[,intersect(gene_XW,colnames(EAS_HGVSp))]
TCGA_HGVSp=as.data.frame(t(read.csv('anno_maf/TCGA_HGVSp.csv',check.names = F)))
colnames(TCGA_HGVSp)=TCGA_HGVSp[1,]
TCGA_HGVSp=TCGA_HGVSp[-1,]
TCGA_HGVSp=TCGA_HGVSp[,intersect(gene_XW,colnames(TCGA_HGVSp))]

drivers=unique(c(colnames(XW_HGVSp),colnames(EAS_HGVSp),colnames(TCGA_HGVSp)))
for (i in setdiff(drivers,colnames(XW_HGVSp))){
  XW_HGVSp[,i]='WT'
}
XW_HGVSp$cohort='XW'
for (i in setdiff(drivers,colnames(EAS_HGVSp))){
  EAS_HGVSp[,i]='WT'
}
EAS_HGVSp$cohort='EAS'
for (i in setdiff(drivers,colnames(TCGA_HGVSp))){
  TCGA_HGVSp[,i]='WT'
}
TCGA_HGVSp$cohort='TCGA'

cohort3_HGVSp=rbind(XW_HGVSp,EAS_HGVSp,TCGA_HGVSp)  
cohort3_HGVSp$value=1
cohort3_HGVSp$cohort=factor(cohort3_HGVSp$cohort,level=c('XW','EAS','TCGA'))
EGFR_SBS4=function(df,cc){
  x=df@data
Y=x[x$Tumor_Sample_Barcode %in% SBS4_immune_matrix$Tumor_Sample_Barcode,]
Y1=Y[Y$Hugo_Symbol=='EGFR',]
Y1$cohort=cc
return(Y1)
}




EGFR_SBS4_XW=EGFR_SBS4(all_mut,'XW')
EGFR_SBS4_EAS=EGFR_SBS4(anno_maf_EAS,'EAS')
EGFR_SBS4_TCGA=EGFR_SBS4(anno_maf_TCGA,'TCGA')

EGFR_SBS4_3cohort=rbind(EGFR_SBS4_XW,EGFR_SBS4_EAS,EGFR_SBS4_TCGA)
EGFR_SBS4_3cohort$value=1
EGFR_SBS4_3cohort$ppos=as.numeric(str_extract(EGFR_SBS4_3cohort$HGVSp_Short, "\\d+"))
EGFR_SBS4_3cohort$structure <- ifelse(EGFR_SBS4_3cohort$ppos >= 718 & EGFR_SBS4_3cohort$ppos <= 726, "p-loop",
                ifelse(EGFR_SBS4_3cohort$ppos < 718, "far_p_loop",
                       ifelse(EGFR_SBS4_3cohort$ppos >= 767 & EGFR_SBS4_3cohort$ppos <= 779, "c-helix",
                              ifelse(ifelse(EGFR_SBS4_3cohort$ppos > 726 & EGFR_SBS4_3cohort$ppos < 767, "inter",
                                            ifelse(EGFR_SBS4_3cohort$ppos > 779,'far-loop',''))))))


cohort3_HGVSp_SBS4=cohort3_HGVSp[intersect(row.names(cohort3_HGVSp),SBS4_immune_matrix$Tumor_Sample_Barcode),]
cohort3_HGVSp_SBS4$typical='Atypical'
cohort3_HGVSp_SBS4$typical[cohort3_HGVSp_SBS4$EGFR %in% c("p.T790M",
                                                          "p.L858R",
                                                          "p.E746_A750del",
                                                          "p.L747_T751del",
                                                          "p.L747_S752del")]='Classical'
cohort3_HGVSp_SBS4$typical[cohort3_HGVSp_SBS4$EGFR=='WT']='WT'
cohort3_HGVSp_SBS4$typical[grep('/',cohort3_HGVSp_SBS4$EGFR)]='Double'
cohort3_HGVSp_SBS4$typical=factor(cohort3_HGVSp_SBS4$typical,levels = c('Classical','Atypical','Double','WT'))
cohort3_HGVSp_SBS4$EGFR=factor(cohort3_HGVSp_SBS4$EGFR,
                               level=c("p.T790M",
                                       "p.L858R",
                                       "p.E746_A750del",
                                       "p.L747_T751del",
                                       "p.L747_S752del",
                                       setdiff(cohort3_HGVSp_SBS4$EGFR,c("p.T790M",
                                                                         "p.L858R",
                                                                         "p.E746_A750del",
                                                                         "p.L747_T751del",
                                                                         "p.L747_S752del",'WT')),'WT'))

p1=ggplot(cohort3_HGVSp_SBS4, aes(x = cohort, y = value, fill = typical)) +  
  #geom_bar(stat = 'identity', position = 'fill', width = 0.75) +  
  geom_bar(position = "stack", stat = "identity",width = 0.75) + 
  scale_fill_manual(values = c("#ec8181","#d2da93",'grey',"#5196d5","#00ceff","#ff630d","#35978b", 
                               "#e5acd7","#77aecd","#dfc6a5","#e50719", 
                               "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
                               "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                               'steelblue','yellow')) +  #
  labs(x='Cohort',y='Frequency')+
  theme_bw()+
  theme(#legend.position='none',
    axis.line = element_line(size=0.25),
    axis.title = element_text(size=18,
                              #face="bold",
                              color="black"),
    panel.grid.major = element_blank(),   # 删除主要网格线
    panel.grid.minor = element_blank(),   # 删除次要网格线
    panel.border = element_rect(color = "black", size = 0.25, fill = NA),
    axis.text = element_text(size=18,color = "black"),
    plot.title=element_text(size=12,
                            face="bold",
                            color="black", #颜色
                            hjust=0.5, #调整位置，正中间
                            lineheight=0.5)
    #plot.background = element_rect(color = "blue", size = 2)
  )
p2=ggplot(cohort3_HGVSp_SBS4, aes(x = cohort, y = value, fill = EGFR)) +  
  geom_bar(stat = 'identity', position = 'fill', width = 0.75) +  
  #geom_bar(position = "stack", stat = "identity",width = 0.75) + 
  scale_fill_manual(values = c("#ec8181","#d2da93",'grey',"#5196d5","#00ceff","#ff630d","#35978b", 
                               "#e5acd7","#77aecd","#dfc6a5","#e50719", 
                               "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
                               "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                               'steelblue','yellow')) +  #
  labs(x='Cohort',y='Frequency')+
  theme_bw()+
  theme(#legend.position='none',
    axis.line = element_line(size=0.25),
    axis.title = element_text(size=18,
                              #face="bold",
                              color="black"),
    panel.grid.major = element_blank(),   # 删除主要网格线
    panel.grid.minor = element_blank(),   # 删除次要网格线
    panel.border = element_rect(color = "black", size = 0.25, fill = NA),
    axis.text = element_text(size=15,color = "black"),
    plot.title=element_text(size=12,
                            face="bold",
                            color="black", #颜色
                            hjust=0.5, #调整位置，正中间
                            lineheight=0.5)
    #plot.background = element_rect(color = "blue", size = 2)
  )
 
p1|p2 
EGFR_SBS4_3cohort$Exon=paste0('Exon ',gsub('/28','',EGFR_SBS4_3cohort$Exon_Number))
EGFR_SBS4_3cohort$Exon[!(EGFR_SBS4_3cohort$Exon %in% c('Exon 18',"Exon 19","Exon 20","Exon 21"))]='Other'
EGFR_SBS4_3cohort$cohort=factor(EGFR_SBS4_3cohort$cohort,levels = c('XW','EAS','TCGA'))
p3=ggplot(EGFR_SBS4_3cohort, aes(x = cohort, y = value, fill = Exon)) +  
  #geom_bar(stat = 'identity', position = 'fill', width = 0.75) +  
  geom_bar(position = "stack", stat = "identity",width = 0.75) + 
  scale_fill_manual(values = c("#8DD3C7","#FFFFB3",'#BEBADA',"#FB7F72",'steelblue',"#00ceff","#ff630d","#35978b", 
                               "#e5acd7","#77aecd","#dfc6a5","#e50719", 
                               "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
                               "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                               'yellow')) +  #
  labs(x='Cohort',y='Frequency')+
  theme_bw()+
  theme(#legend.position='none',
    axis.line = element_line(size=0.25),
    axis.title = element_text(size=18,
                              #face="bold",
                              color="black"),
    panel.grid.major = element_blank(),   # 删除主要网格线
    panel.grid.minor = element_blank(),   # 删除次要网格线
    panel.border = element_rect(color = "black", size = 0.25, fill = NA),
    axis.text = element_text(size=18,color = "black"),
    plot.title=element_text(size=12,
                            face="bold",
                            color="black", #颜色
                            hjust=0.5, #调整位置，正中间
                            lineheight=0.5)
    #plot.background = element_rect(color = "blue", size = 2)
  )
p1|p3
###纠正TCGA筛选
X1=read.csv('tcga.csv')
x=X1[X1$submitter_id.samples %in% cohort4_mixcopdata_Nofemale$Tumor_Sample_Barcode,]
x=x[is.na(x$stopped_smoking_year)==F,]
x$stop_years=x$year_of_initial_pathologic_diagnosis-x$stopped_smoking_year
x1=x[x$stop_years>=10&x$gender.demographic=='female',]
x2=immune_matrix[immune_matrix$Tumor_Sample_Barcode %in% x1$submitter_id.samples,]


##signature_gsva------------------------------
signature_gsva$Clock_like_signature_gsva=signature_gsva$SBS1+signature_gsva$SBS5

plot_data=signature_gsva[,c("cohort",colnames(signature_gsva)[grep('_signature',colnames(signature_gsva))])]
X=data.frame(Cohort=factor(plot_data[,1],levels = c('XW','EAS','TCGA')),
             row.names = row.names(plot_data))

# 构建颜色映射
ann_ <- list(
  Cohort = list(
    Cohort = c("EAS" = "#E3C06C", "XW" = "#D07358", "TCGA" = "#2A978A") # 替换为实际分类和颜色
  )
)

#pheatmap(t(plot_data[,c(2:11,13,15,16)]),
pheatmap(t(plot_data[,c(2,6,8:11,16)]),
         scale = 'rows',
         color = colorRampPalette(colors = c("grey", "#d56e5e"))(100),
         cluster_cols = T,
         cluster_rows = T,
         #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
         clustering_distance_rows = "canberra",
         #"ward.D", "ward.D2", "single", "complete", 
         #"average" (= UPGMA), "mcquitty" (= WPGMA), 
         #"median" (= WPGMC) or "centroid" (= UPGMC)
        clustering_method = "complete",
         show_colnames = F,
         
         
         annotation_col =X,
         annotation_colors = ann_$Cohort 
)
