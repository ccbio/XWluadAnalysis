
#/bin/bash
#PBS -N CNV_call_4
#PBS -q cypQueue
#PBS -l nodes=node25:ppn=8
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
out_path=${project}AligementOut/WES-4/
work_path=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/
#SampleList=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch1_group_WES_pairwise.csv
tmp_dir=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/step

tumorlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch4_tumor_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $tumorlist))'
sample_tumor=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

normallist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch4_normal_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $normallist))'
sample_normal=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

IDlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch4_ID_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $IDlist))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

normal=${out_path}${sample_normal}/${sample_normal}.sorted.markdup.BQSR_hg38_bed.bam
tumor=${out_path}${sample_tumor}/${sample_tumor}.sorted.markdup.BQSR_hg38_bed.bam

source activate sciclone_env

mkdir ${work_path}/cnv_sequenza
outdir=${work_path}/cnv_sequenza
mkdir ${outdir}/$sampleID

sequenza-utils bam2seqz -n ${normal} \
	-t ${tumor}  \
	--fasta ${genome} \
	-gc ${gc} \
	-C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
	--parallel 12 \
	-o ${outdir}/$sampleID/${sampleID}.seqz.gz

cd $outdir/$sampleID
source activate gatk
for i in {1..22}
do
sequenza-utils seqz_binning --seqz ${sampleID}_chr${i}.seqz.gz -w 50 -o ${sampleID}_chr${i}.small.seqz.gz
done
sequenza-utils seqz_binning --seqz ${sampleID}_chrX.seqz.gz -w 50 -o ${sampleID}_chrX.small.seqz.gz
sequenza-utils seqz_binning --seqz ${sampleID}_chrY.seqz.gz -w 50 -o ${sampleID}_chrY.small.seqz.gz

$(> temp)
for i in `ls ${sampleID}_chr*.small.seqz.gz`;do
	echo -n "$i " >> temp
done
h2=$(cat temp)
h2=${h2%*,}
#rm ${sampleID}.seqz.gz
#rm ${sampleID}.seqz.gz.tbi
zcat $h2 | gawk '{if (NR!=1 && $1 != "chromosome"){print $0}}'| bgzip > ${sampleID}_small.seqz.gz
tabix -f -s 1 -b 2 -e 2 -S 1 ${sampleID}_small.seqz.gz


source activate r4_py37_env
Rscript /data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/workflow/CNV/sequenza.R ${sampleID}_small.seqz.gz ${sampleID} ./  
