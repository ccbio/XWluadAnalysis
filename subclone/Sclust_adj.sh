#/bin/bash
#PBS -N Sclust 
#PBS -q cypQueue
#PBS -l nodes=node11:ppn=2
#PBS -l mem=30gb
#PBS -o QsubLog/
#PBS -e QsubLog/

project=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/
ref_path=/data_group/cunyupeng/hg38/hg38.fa
out_path=${project}AligementOut/

#tumorlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/workflow/Sclust/1.5_0.08_3.5_100_1e-7_PASS_1.2_tumor_up.txt
#IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $tumorlist))'
#sample_tumor=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sample_tumor=1224566-TT7_FDHE19H000620-1a-A40-A61 #`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#tumor=1196553-TT6_FDHE19H000589-1a-A15-A51
#normallist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/workflow/Sclust/1.5_0.08_3.5_100_1e-7_PASS_1.2_normal_up.txt
#IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $normallist))'
#sample_normal=1224566-TP7_FDHE19H000626-1a-A39-A51 #`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sample_normal=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#normal=1196553-TP6_FDHE19H000592-1a-A14-A30

IDlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/workflow/Sclust/1.5_0.08_3.5_100_1e-7_PASS_0.85_ID_up.txt
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $IDlist))'
#sampleID=1224566 #`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sampleID=1196553

#sample_normal=WES-2/1215795-TP3_FDHE19H000625-1a-A21-A4/1215795-TP3_FDHE19H000625-1a-A21-A4
#sample_tumor=WES-2/1215795-TT3_FDHE19H000625-1a-A22-A36/1215795-TT3_FDHE19H000625-1a-A22-A36
#sampleID=1249132

normal=${out_path}${sample_normal}.sorted.markdup.BQSR_hg38_bed.bam
tumor=${out_path}${sample_tumor}.sorted.markdup.BQSR_hg38_bed.bam
echo ${tumor}

source activete base
cd ${project}

#if [ -d "/Sclust/" ];then
#  rm -r  Sclust
#else
#  mkdir Sclust
#fi

#mkdir ${project}/Sclust/${sampleID}/
Sclust=/opt/service/Sclust/bin/Sclust
#Sclust=/data_group/cunyupeng/sws/Sclust/bin/Sclust
workp=${project}/somatic/Sclust/$sampleID/
cd $workp

#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr1
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr2
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr3
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr4
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr5
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr6
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr7
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr8
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr9
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr10
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr11
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr12
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr13
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr14
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr15
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr16
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr17
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr18
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr19
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr20
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr21
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chr22
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chrX
#${Sclust} bamprocess -t ${tumor} -n ${normal} -o ${workp}$sampleID -part 1 -build hg38 -r chrY
#
#Sclust bamprocess -i ${workp}${sampleID} -o ${workp}${sampleID}
#
#$Sclust cn -rc ${sampleID}_rcount.txt -snp ${sampleID}_snps.txt -vcf ${sampleID}_mutect2.vcf -o ${sampleID}_mutect2
#$Sclust cn -rc ${sampleID}_rcount.txt -snp ${sampleID}_snps.txt -vcf ${sampleID}_varscan.vcf -o ${sampleID}_varscan
#$Sclust cn -rc ${sampleID}_rcount.txt -snp ${sampleID}_snps.txt -vcf ${sampleID}_strelka.vcf -o ${sampleID}_strelka
#rm all.vcf
#for i in *vcf;do tail -n  +2 $i >>all.vcf; done
#rm *chr*	
cd ${sampleID}
minp=1.5
minpu=0.01
maxp=3.5
nbins=100
lambda=1e-7

$Sclust cn -rc ${sampleID}_rcount.txt -snp ${sampleID}_snps.txt -vcf all_select__pass_2.vcf \
	-minp ${minp} -minpu ${minpu} -maxp ${maxp} \
	-o ${sampleID}_${minp}_${minpu}_${maxp}_${nbins}_${lambda}_PASS 

$Sclust cluster -i ${sampleID}_${minp}_${minpu}_${maxp}_${nbins}_${lambda}_PASS \
	-nbins ${nbins} -lambda ${lambda}

