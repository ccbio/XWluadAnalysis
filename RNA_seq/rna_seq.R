library(sva)
gene_idint=function(df,gene_int){
  row.names(df)=gene_int
  df1=df[,-1]
  return(df1)
}

path='xxxxx'
tcga_lusc=read.table(paste0(path,'TCGA-LUSC.htseq_counts.tsv'),header = 1,check.names = F)

path='xxxxx'
tcga_luad=read.table(paste0(path,'TCGA-LUAD.htseq_counts.tsv'),header = 1,check.names = F)   

EAS_luad=read.table(header = 1,'GSK_RSEM_rerun_expCounts_172_Tumor.tsv')
XW_luad=read.table(header = 1,'XW_RSEM_expCounts_150.tsv')
EAS_luad_normal=read.table(header = 1,'GSK_RSEM_rerun_expCounts_88_Normal.tsv')


gene_trans=EAS_luad[,c(1:3)]
gene_trans$Ensembl.ID=gsub("\\.[0-9]*",'',gene_trans$Ensembl.ID)


EAS_gene=gsub("\\.[0-9]*",'',EAS_luad$Ensembl.ID)
EAS_luad1=gene_idint(EAS_luad,EAS_gene)
EAS_luad1=EAS_luad1[,-c(1,2)]
EAS_gene_normal=gsub("\\.[0-9]*",'',EAS_luad_normal$Ensembl.ID)
EAS_luad_normal1=gene_idint(EAS_luad_normal,EAS_gene)
EAS_luad_normal1=EAS_luad_normal1[,-c(1,2)]



tcga_luad_gene=gsub("\\.[0-9]*",'',tcga_luad$Ensembl_ID)
tcga_luad1=gene_idint(tcga_luad,tcga_luad_gene)
tcga_lusc_gene=gsub("\\.[0-9]*",'',tcga_lusc$Ensembl_ID)
tcga_lusc1=gene_idint(tcga_lusc,tcga_lusc_gene)
XW_gene=gsub("\\.[0-9]*",'',XW_luad$gene_id)
XW_luad1=gene_idint(XW_luad,XW_gene)
GENE_list=list(EAS_gene=EAS_gene,
               tcga_luad_gene=tcga_luad_gene,
               tcga_lusc_gene=tcga_lusc_gene,
               XW_gene=XW_gene)

gene_id=intersect(intersect(intersect(GENE_list$EAS_gene,GENE_list$tcga_luad_gene),GENE_list$XW_gene),GENE_list$tcga_lusc_gene)

XW_luad2=XW_luad1[gene_id,]

EAS_luad2=EAS_luad1[gene_id,]
EAS_luad_normal2=EAS_luad_normal1[gene_id,]
tcga_luad2=tcga_luad1[gene_id,]
tcga_luad2=tcga_luad2
tcga_lusc2=tcga_lusc1[gene_id,]
tcga_lusc2=tcga_lusc2
cohort3_genecount=cbind(log2(XW_luad2+1),log2(EAS_luad2+1),tcga_luad2)


group_tcga_luad=data.frame(sample=colnames(tcga_luad2),cohort='TCGA_LUAD')
#group_tcga_luad=group_tcga_luad[-1,]
group_tcga_lusc=data.frame(sample=colnames(tcga_lusc2),cohort='TCGA_LUSC')
#group_tcga_lusc=group_tcga_lusc[-1,]
group_XW_luad=data.frame(sample=colnames(XW_luad2),cohort='XW')
#group_XW_luad=group_XW_luad[-1,]
group_EAS_luad=data.frame(sample=colnames(EAS_luad2),cohort='EAS')
#group_EAS_luad=group_EAS_luad[-1,]

cohort3_GROUP=rbind(group_XW_luad,group_EAS_luad,group_tcga_luad)


tumor_sample=function(df){
  NoFemale_XW<-readRDS(df)
  NoFemale_XW=as.data.frame(NoFemale_XW)
  xx=row.names(NoFemale_XW[NoFemale_XW$ V1==1,])
  return(xx)
}
NoFemale_XW<-tumor_sample("matrix_XW.rdata")
NoFemale_EAS<-tumor_sample("matrix_EAS.rdata")
NoFemale_TCGA<-tumor_sample("matrix_LUAD.rdata")


cohort3_genecount_combat=ComBat(cohort3_genecount,batch = cohort3_GROUP$cohort)
cohort3_genecount_nofemale=cohort3_genecount_combat[,c(NoFemale_XW,NoFemale_EAS,NoFemale_TCGA)]
cohort3_genecount_nofemale=as.data.frame(cohort3_genecount_nofemale)

x=gene_trans[gene_trans$Gene.type=='protein-coding',]
x=na.omit(x)
row.names(x)=x$Ensembl.ID
y=x[intersect(x$Ensembl.ID,row.names(cohort3_genecount_nofemale)),]

cohort3_genecount_nofemale3=cohort3_genecount_nofemale[intersect(x$Ensembl.ID,row.names(cohort3_genecount_nofemale)),]
cohort3_genecount_nofemale3$gene_id=y[,2]
cohort3_genecount_nofemale4=aggregate(.~gene_id,data=cohort3_genecount_nofemale3,sum)
row.names(cohort3_genecount_nofemale4)=cohort3_genecount_nofemale4[,1]
cohort3_genecount_nofemale4=cohort3_genecount_nofemale4[,-1]

all_stsvm_gene=union(union(GENE_list$TCGA,GENE_list$EAS),GENE_list$XW)
cohort3_genecount_nofemale_stsvm=cohort3_genecount_nofemale4[all_stsvm_gene,]

saveRDS(cohort3_genecount_nofemale_stsvm,file = 'cohort3_genecount_nofemale_stsvm.RDS')

####SNF

##XW
x1=str_split_fixed(NoFemale_XW, "_", n = 3)[,2]
x1=str_split_fixed(x1, "\\.", n = 2)[,1]
x1=gsub('s','',x1)

tumor_sample_matrix=function(df){
  NoFemale_XW<-readRDS(df)
  NoFemale_XW=as.data.frame(NoFemale_XW)
  xx=row.names(NoFemale_XW[NoFemale_XW$ V1==1,])
  return(NoFemale_XW[xx,-1])
}
NoFemale_XW_EXPR<-tumor_sample_matrix("matrix_XW.rdata")
NoFemale_XW_CNV<-read.csv(sep = '\t','cnv_sequenza.mat',check.names = F)
XW_mut <- read.maf('SNP_hg38_all_mut.maf')
NoFemale_XW_SNV<-mutCountMatrix(XW_mut)

x2=intersect(intersect(colnames(NoFemale_XW_SNV),x1),colnames(NoFemale_XW_CNV))
XW_list_SNF=list(RNA=NoFemale_EAS_EXPR,CNV=NoFemale_XW_CNV,SNV=as.data.frame(EAS_mutconut))

##EAS

NoFemale_EAS_EXPR<-cohort3_genecount_nofemale4[,NoFemale_EAS]
NoFemale_EAS_CNV<-read.csv(sep = '\t','WES/subtype/subtype/XW/movics/EAS/EAS_all.mat',check.names = F)
x1=intersect(NoFemale_EAS,colnames(NoFemale_EAS_CNV))
NoFemale_EAS_EXPR=NoFemale_EAS_EXPR[,x1]
NoFemale_XW_CNV=NoFemale_XW_CNV[,x1]
zww_mut <- read.csv(sep='\t','snv_indel.maf')
EAS_maf=read.maf(zww_mut)
EAS_mutconut=mutCountMatrix(EAS_maf)
EAS_mutconut=EAS_mutconut[,x1]
EAS_list_SNF=list(RNA=NoFemale_EAS_EXPR,CNV=NoFemale_XW_CNV,SNV=as.data.frame(EAS_mutconut))
pdatzww=read.csv('GIS031.clinical.patient.txt',sep='\t')
mixcopdata_ZWW=pdatzww[,c("Tumor_Sample_Barcode",
                          "COHORT","AGE",
                          "GENDER","SMOKING_STATUS",
                          "STAGE")]
A=mixcopdata_ZWW[mixcopdata_ZWW$GENDER=='Female'&mixcopdata_ZWW$SMOKING_STATUS=='No',]
a=intersect(A$Tumor_Sample_Barcode,x1)


TCGA_mut=read.csv(sep = '\t','TCGA-LUAD.mutect2_snv.new.tsv',header = T)