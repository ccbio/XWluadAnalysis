#!/bin/bash
#PBS -N Rsem
#PBS -o QsubLog/
#PBS -e QsubLog/
#PBS -q cypQueue
#PBS -l mem=5gb
#PBS -l nodes=node23:ppn=4

################################
### Variable information ###
MainFile=Hg38RNAseq

### Project information ###
ProjectDir=/data_group/cunyupeng/huling/daiyi/Project/${MainFile}
SampleList=${ProjectDir}/Meta/remaind.idx
RawReads=${ProjectDir}/Raw

PBSLog=${ProjectDir}/Workflow/QsubLog

CleanReads=${ProjectDir}/CleanReads
FastQC=${ProjectDir}/FastQCDir
FastpDir=${ProjectDir}/FastpDir

################################
### Main parameters ###
ReferenceGenome=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/ref/GCF_000001405.31_GRCh38.p5_genomic.primary_assembly.fna
ReferenceGtf=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/ref/GCF_000001405.31_GRCh38.p5_genomic.gtf
RequiredCPU=4


################################
### Tools ###
Fastp=fastp
Fastqc=fastqc

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sample_idx=${sample_idx}
sampleID=`echo ${ARRAY[((${sample_idx}-1))]}`

#sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

################################
### Main process ###

cd ${ProjectDir}/Star
source activate rex_env

###Build RSEM references using RefSeq###
#${Rsem} ${ReferenceGenome} ${ProjectDir}/Homo_sapiens_assembly38.primary_assembly.fa

rsem-calculate-expression --paired-end --star --star-gzipped-read-file -p 4 ${CleanReads}/${sampleID}_1.fq.gz ${CleanReads}/${sampleID}_2.fq.gz ${ProjectDir}/Index/hg38_star ${sampleID}
