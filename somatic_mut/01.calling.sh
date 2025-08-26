#/bin/bash
#PBS -N mapping 
#PBS -q cypQueue
#PBS -l nodes=1:ppn=8
#PBS -l mem=30gb
#PBS -o QsubLog/
#PBS -e QsubLog/
project=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/
ref_path=/data_group/cunyupeng/hg38/hg38.fa
seq_path=${project}metadata/WES/
out_path=${project}AligementOut/
exome_path=/data_group/cunyupeng/hg38/hg38_exon.bed
dbsnp1kG=/data_group/cunyupeng/hg38/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz
dbindel1kG=/data_group/cunyupeng/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
dbindelhg38=/data_group/cunyupeng/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
dbsnp138=/data_group/cunyupeng/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz
afhg38=/data_group/cunyupeng/hg38/somatic-hg38_af-only-gnomad.hg38.vcf.gz
hapmaphg38=/data_group/cunyupeng/hg38/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz
SampleList=${project}workflow/samplelist.txt

IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sample=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sample=69

q1=${seq_path}${sample}_1.clean.fq.gz
q2=${seq_path}${sample}_2.clean.fq.gz
mkdir $out_path/$sample
out_file=$out_path/$sample

## clean
source activate reseq_env
fastp -i $q1 -I $q2 -o ${out_file}/${sample}_R1.clean.fastq.gz -O ${out_file}/${sample}_R2.clean.fastq.gz



bwa mem -t 24 -M -R "@RG\tID:${sample}\tPL:illumina\tLB:WES\tSM:${sample}" ${ref_path} ${out_file}/${sample}_R1.clean.fastq.gz ${out_file}/${sample}_R2.clean.fastq.gz > ${out_file}/${sample}.sam

samtools view -Sb ${out_file}/${sample}.sam > ${out_file}/${sample}.bam

samtools sort -@ 24 -O bam -o ${out_file}/${sample}.sorted.bam ${out_file}/${sample}.bam

gatk MarkDuplicates -I ${out_file}/${sample}.sorted.bam -O ${out_file}/${sample}.sorted.markdup.bam -M ${out_file}/${sample}.sorted.markdup.txt

samtools index ${out_file}/${sample}.sorted.markdup.bam

gatk CreateSequenceDictionary -R ${ref_path} -O ${out_file}/hg38.dict

samtools faidx ${ref_path}

gatk BaseRecalibrator -R ${ref_path} -I ${out_file}/${sample}.sorted.markdup.bam \
                                     -L ${exome_path} \
                                     --use-original-qualities \
                                     -O ${out_file}/${sample}data.table \
                                    --known-sites $dbsnp1kG \
                                    --known-sites $dbindel1kG \
                                    --known-sites $dbsnp138

gatk ApplyBQSR -R ${ref_path} \
               -I ${out_file}/${sample}.sorted.markdup.bam \
               -L ${exome_path} \
               --bqsr-recal-file ${out_file}/${sample}data.table \
               -O ${out_file}/${sample}.sorted.markdup.BQSR.bam

#germline:
gatk HaplotypeCaller -R ${ref_path} \
                     -L ${exome_path} \
                     -I ${out_file}/${sample}.sorted.markdup.BQSR.bam \
                     -O ${out_file}/${sample}.vcf
