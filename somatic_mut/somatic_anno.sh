#/bin/bash
#PBS -N somatic_anno 
#PBS -q cypQueue
#PBS -l nodes=node09:ppn=4
#PBS -l mem=30gb
#PBS -o QsubLog/somatic
#PBS -e QsubLog/somatic

work_path=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/somatic

samplelist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/all_ID_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $samplelist))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sampleID=1196553 #`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
outfile=${work_path}/${sampleID}
samplevcf=${outfile}/${sampleID}_filtered.vcf.gz \
#samplevcf=${outfile}/strelka_pon1/results/variants/somatic.snvs.vcf.gz
source activate reseq_env
gatk Funcotator --variant ${samplevcf} \
                --reference  /data_group/cunyupeng/hg38/hg38.fa\
                --ref-version hg38 \
                --data-sources-path /data_group/cunyupeng/public_data/funcotator_dataSources.v1.7.20200521s  \
                --output ${outfile}/${sampleID}_mutect_funcotator.maf \
                --output-file-format MAF
