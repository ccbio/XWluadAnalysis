#/bin/bash
#PBS -N CNV_ascat 
#PBS -q cypQueue
#PBS -l nodes=node12:ppn=4
#PBS -l mem=30gb
#PBS -o QsubLog/
#PBS -e QsubLog/

project=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/
batch=4
out_path=${project}AligementOut/WES-${batch}/

tumorlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch${batch}_tumor_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $tumorlist))'
sample_tumor=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

normallist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch${batch}_normal_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $normallist))'
sample_normal=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

IDlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch${batch}_ID_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $IDlist))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

normal=${out_path}${sample_normal}/${sample_normal}.sorted.markdup.BQSR_hg38_bed.bam
tumor=${out_path}${sample_tumor}/${sample_tumor}.sorted.markdup.BQSR_hg38_bed.bam


source activate r4.2_env
export PATH=$PATH:/data_group/cunyupeng/sunyh_home/resource/alleleCount/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data_group/cunyupeng/sunyh_home/resource/alleleCount/lib

pathw=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/

ascatpath=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/cnv_ascat/

echo $tumor 
echo $normal 
echo $sample_tumor 
echo $sample_normal 
echo $sampleID
 
mkdir $ascatpath/$sampleID/
cd $ascatpath/$sampleID/
echo $ascatpath/$sampleID/
Rscript $pathw/workflow/CNV/ASCAT.R $tumor $normal $sample_tumor $sample_normal ./ $sampleID
