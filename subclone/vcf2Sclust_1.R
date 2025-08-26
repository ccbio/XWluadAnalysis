mu2sclust <- function(ff= "k36.out",outdir='output',samplepre='')
{	
  mvcf=read.csv(ff,sep="\t",skip=523,header=T)
  vv=data.frame()
  for(i in 1:nrow(mvcf))
  {
    tmp=mvcf[i,]
    tt=strsplit(as.character(tmp[,10]),":")[[1]]
    tn=strsplit(as.character(tmp[,11]),":")[[1]]	
    ti=strsplit(as.character(tmp[,8]),";")[[1]]
    #dp=ti[2]
    dp=paste("DP",tt[4],sep="=") #paste("DP",strsplit(tt[2],",")[[1]][4],sep="=")
    af= paste("AF",tt[3],sep="=")
    #dpn=paste("DP_N", strsplit(tn[2],",")[[1]][2],sep="=") 
    dpn=paste("DP_N",tn[4],sep="=") 
    afn= paste("AF_N",tn[3],sep="=")
    infor= paste(paste(dp,paste(af,dpn,sep=";"),sep=";"), afn,sep=";")
    vv=rbind(vv,data.frame(INFO=infor))
  }
  vvv= cbind(mvcf[,1:7],INFO=vv)
  write.table(cbind(mvcf[,1:7],INFO=vv),file=paste(outdir,samplepre,"_mutect2.vcf",sep=""),sep="\t",quote = F,row.names=F)
  return(vvv)
}

vs2sclust <- function(ff= "k36.out")
{ 
  mvcf=read.csv(ff,sep="\t",skip=17,header=T)
  vv=data.frame()
  for(i in 1:nrow(mvcf))
  {
    tmp=mvcf[i,]
    tn=strsplit(as.character(tmp[,10]),":")[[1]]
    tt=strsplit(as.character(tmp[,11]),":")[[1]]	
    ti=strsplit(as.character(tmp[,8]),";")[[1]]
    #dp=ti[2]
    dp=paste("DP",tt[3],sep="=") #paste("DP",strsplit(tt[2],",")[[1]][4],sep="=")
	vaf=as.numeric (sub ("%","",tt[6]))/100
	
    af= paste("AF",vaf,sep="=")
    #dpn=paste("DP_N", strsplit(tn[2],",")[[1]][2],sep="=") 
    dpn=paste("DP_N",tn[3],sep="=") 
    vafn=as.numeric (sub ("%","",tn[6]))/100
	afn= paste("AF_N",vafn,sep="=")
    infor= paste(paste(dp,paste(af,dpn,sep=";"),sep=";"), afn,sep=";")
    vv=rbind(vv,data.frame(INFO=infor))
  }
  ss=strsplit(ff,".f")[[1]][1]
  
  vvv= cbind(mvcf[,1:7],INFO=vv)
  #write.table(cbind(mvcf[,1:7],INFO=vv),file=paste(ss,"_varscan.vcf",sep=""),sep="\t",quote = F,row.names=F)
  return(vvv)
}

## strelka2sclsutVCF
##
skvs2sclust <- function(mvcf= mvcf)
{
  nn=c("A","C","G","T")
  vv=data.frame()
  vv1=data.frame()
  for(i in 1:nrow(mvcf))
  {if(i%%10000==1)#print(i)
    tmp=mvcf[i,]
    #get read count for each NN
    #
    ni=strsplit(as.character(tmp$NORMAL),":")[[1]][5:8]; ti=strsplit(as.character(tmp$TUMOR),":")[[1]][5:8]
    nii=rep(0,4); for(i in 1:length(ni))nii[i]=as.numeric(strsplit(ni[i],",")[[1]][1])
    tii=rep(0,4); for(i in 1:length(ti))tii[i]=as.numeric(strsplit(ti[i],",")[[1]][1])
    ##get read count for Ref/ALT in normal/tumor
    ref=strsplit(as.character(tmp$REF),"")[[1]]; alt=strsplit(as.character(tmp$ALT),"")[[1]]
    tumor_rRef=sum(tii[which(nn %in% ref)])
    tumor_rAlt=sum(tii[which(nn %in% alt)])
    normal_rRef=sum(nii[which(nn %in% ref)])
    normal_rAlt=sum(nii[which(nn %in% alt)])
    
    #allele frequency of tumor/normal
    dpt1=tumor_rRef+tumor_rAlt
    dpn1=normal_rRef+normal_rAlt
    
    dp=paste("DP",dpt1,sep="=") 
    af= paste("AF",tumor_rAlt/dpt1,sep="=")
    dpn=paste("DP_N",dpn1,sep="=") 
    afn= paste("AF_N",normal_rAlt/dpn1,sep="=")
    
    infor= paste(paste(dp,paste(af,dpn,sep=";"),sep=";"), afn,sep=";")
    
    vv1=rbind(vv1,data.frame(af_t=tumor_rAlt/dpt1,af_n=normal_rAlt/dpn1,dp_t=dpt1,dp_n=dpn1))
    vv=rbind(vv,data.frame(INFO=infor))
    
  }
  vvv= cbind(mvcf[,1:7],INFO=vv)
  vvv1= cbind(mvcf[,1:7],vv1)
  
  chrs=paste("chr",c(1:22,"X","Y"),sep="")
  vvv=vvv[which(vvv$X.CHROM %in% chrs),]
  vv1=vv1[which(vvv$X.CHROM %in% chrs),]
  vvvs=data.frame()
  vvvs1=data.frame()
  for(c in chrs)
  {
    tmp= vvv[which(vvv[,1]==c),]
    tmp=tmp[order(tmp[,2],decreasing = F),]
    vvvs=rbind(vvvs,tmp)
    tmp1= vvv1[which(vvv1[,1]==c),]
    tmp1= tmp1[order(tmp[,2],decreasing = F),]
    vvvs1=rbind(vvvs1,tmp1)
  }
  
  return(res=list(vcf=vvvs,vcfa=vvvs1))
} 
pairedvcf <- function(vcfsnp="df1",vcfindel="df2",soft='varscan',outdir='',samplepre='' ){
  
  if (soft=='strelka'){
    v1=read.csv(vcfsnp,sep="\t",header=T,skip=489)
	print(colnames(v1))    
	v2=read.csv(vcfindel,skip=490,sep="\t",header=T)
    v2f=skvs2sclust(mvcf=v2)
    v1f=skvs2sclust(mvcf=v1)
    vv=rbind(v1f$vcf,v2f$vcf)
    tmp=paste(outdir,samplepre,"_strelka.vcf",sep="")
    write.table(vv,file=tmp,sep="\t",quote = F,row.names=F)
  } else if (soft=='varscan'){
    v2f=vs2sclust(ff=vcfsnp)
    v1f=vs2sclust(ff=vcfindel)
    vv=rbind(v1f,v2f)
    tmp=paste(outdir,samplepre,"_varscan.vcf",sep="")
    write.table(vv,file=tmp,sep="\t",quote = F,row.names=F)
  }
}

args <- commandArgs()
input_path=args[6]
output_path=args[7]
sampleID=args[8]
print(args[8])

mut1=paste0(input_path,sampleID,"/",sampleID,'_filtered.vcf')

var1=paste0(input_path,sampleID,"/",sampleID,".snp.Somatic.hc.filter.vcf")
var2=paste0(input_path,sampleID,"/",sampleID,".indel.Somatic.hc.vcf")

sk1=paste0(input_path,sampleID,"/strelka/results/variants/somatic.snvs.vcf")
sk2=paste0(input_path,sampleID,"/strelka/results/variants/somatic.indels.vcf")

#=========try++++++++++
mu2sclust(ff=mut1,outdir=output_path,samplepre=sampleID )

pairedvcf(vcfsnp = var1,vcfindel = var2,soft='varscan',outdir=output_path,samplepre=sampleID )
print('vaescan')
print(sk1)
print(sk2)
pairedvcf(vcfsnp = sk1,vcfindel = sk2,soft='strelka',outdir=output_path,samplepre=sampleID )





