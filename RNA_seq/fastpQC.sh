#!/bin/bash
#PBS -N Rsem_QC
#PBS -o QsubLog/
#PBS -e QsubLog/
#PBS -q cypQueue
#PBS -l mem=5gb
#PBS -l nodes=node12:ppn=4


################################
### Variable information ###
MainFile=Hg38RNAseq

### Project information ###
ProjectDir=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/${MainFile}
SampleList=${ProjectDir}/Meta/batch4.idx
RawReads=${ProjectDir}/Raw

Index=${ProjectDir}/Index
PBSLog=${ProjectDir}/Workflow/QsubLog/RsemBowtie

CleanReads=${ProjectDir}/CleanReads
FastQC=${ProjectDir}/FastQCDir
FastpDir=${ProjectDir}/FastpDir

Bowtie=${ProjectDir}/Bowtie

################################
### Tools ###
Rsem=rsem-prepare-reference
Bowtie=bowtie2

################################
### Main parameters ###
ReferenceGenome=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/ref/GCF_000001405.31_GRCh38.p5_genomic.fna.gz
ReferenceGtf=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/ref/GCF_000001405.31_GRCh38.p5_genomic.gtf
RequiredCPU=4


################################
### Tools ###
Fastp=fastp
Fastqc=fastqc

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
#sample_idx=${sample_idx}
#sampleID=`echo ${ARRAY[((${sample_idx}-1))]}`

sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

################################
### Main process ###
cd ${ProjectDir}

### FastQC for raw reads ###
${Fastqc} -o ${FastQC}/RawFastQC -t ${RequiredCPU} ${RawReads}/${sampleID}_1.clean.fq.gz ${RawReads}/${sampleID}_2.clean.fq.gz

### Quality control ###
${Fastp} -i ${RawReads}/${sampleID}_1.clean.fq.gz -I ${RawReads}/${sampleID}_2.clean.fq.gz -o ${CleanReads}/${sampleID}_1.fq.gz \
        -O ${CleanReads}/${sampleID}_2.fq.gz -W 5 -M 20 -5 -3 -l 50 -w ${RequiredCPU} \
        -j ${FastpDir}/${sampleID}.json -h ${FastpDir}/${sampleID}.html > ${FastpDir}/FastpLog/${sampleID}.qc.log 2>&1

### Fastqc for clean reads ###
${Fastqc} -o ${FastQC}/CleanFastQC -t ${RequiredCPU} ${CleanReads}/${sampleID}_1.fq.gz ${CleanReads}/${sampleID}_2.fq.gz


################################
### adding a finished tag
echo fastp_QC of ${sampleID} is finished at `date`

#qsub -l nodes=node23 node24 node25 -v sample_idx=${i} fastpQC.sh
