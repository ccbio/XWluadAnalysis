

#/bin/bash
#PBS -N CNV_call 
#PBS -q cypQueue
#PBS -l nodes=node06:ppn=4
#PBS -l mem=30gb
#PBS -o QsubLog
#PBS -e QsubLog

gc=/data_group/cunyupeng/sunyh_home/resource/hg38/hg38.gc50Base.wig.gz
bed=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/ref/hg38_exon_chr.bed
#pon=/data_group/cunyupeng/hg38/somatic-hg38_1000g_pon.hg38.vcf.gz
pon=/data_group/cunyupeng/sunyh_home/resource/try/hg38/1000g_pon.hg38.vcf.gz
ge=/data_group/cunyupeng/hg38/somatic-hg38_af-only-gnomad.hg38.vcf.gz
project=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/
genome=/data_group/cunyupeng/hg38/hg38.fa
out_path=${project}AligementOut/WES-1/
work_path=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/
#SampleList=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch1_group_WES_pairwise.csv
tmp_dir=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/step

tumorlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch1_tumor_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $tumorlist))'
sample_tumor=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

normallist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch1_normal_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $normallist))'
sample_normal=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

IDlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch1_ID_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $IDlist))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

normal=${out_path}${sample_normal}/${sample_normal}.sorted.markdup.BQSR_hg38_bed.bam
tumor=${out_path}${sample_tumor}/${sample_tumor}.sorted.markdup.BQSR_hg38_bed.bam

mkdir ${work_path}/cnv_varscan
outdir=${work_path}/cnv_varscan
cd $outdir
mkdir $sampleID
cd $sampleID 
source activate gatk

samtools mpileup -f $genome $normal > ${sampleID}_normal.pileup
samtools mpileup -f $genome $tumor  > ${sampleID}_tumor.pileup

varscan copynumber  ${sampleID}_normal.pileup  ${sampleID}_tumor.pileup ${sampleID}_out
varscan copyCaller ${sampleID}_out.copynumber --output-file ${sampleID}_out.copynumber.called

path=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/workflow/CNV
source activate r4_py37_env
Rscript $path/DNAcopy.R ${sampleID}_out.copynumber.called ${sampleID}_out.copynumber.called.seg

perl $path/mergeSegments.pl ${sampleID}_out.copynumber.called.seg  --ref-arm-sizes $path/inputrefarmsize.txt --output ${sampleID}_out
