#/bin/bash
#PBS -N somatic_call 
#PBS -q cypQueue
#PBS -l nodes=node25:ppn=8
#PBS -l mem=30gb
#PBS -o QsubLog/somatic
#PBS -e QsubLog/somatic

#bed=/data_group/cunyupeng/hg38test/hg38_Exome.bed
bed=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/ref/hg38_Exome.bed
#bed=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/ref/hg38_exon_chr.bed
#pon=/data_group/cunyupeng/hg38/somatic-hg38_1000g_pon.hg38.vcf.gz
pon=/data_group/cunyupeng/sunyh_home/resource/try/hg38/1000g_pon.hg38.vcf.gz
ge=/data_group/cunyupeng/hg38/somatic-hg38_af-only-gnomad.hg38.vcf.gz
project=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/
genome=/data_group/cunyupeng/hg38/hg38.fa
out_path=${project}AligementOut/WES-4/
work_path=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/somatic
#SampleList=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch1_group_WES_pairwise.csv
tmp_dir=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/step

tumorlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch4_tumor_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $tumorlist))'
sample_tumor=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sample_tumor=4824566-TT7_FDHE19H000620-1a-A40-A61 #`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sample_tumor=1196553-TT6_FDHE19H000589-1a-A15-A51
normallist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch4_normal_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $normallist))'
#sample_normal=4824566-TP7_FDHE19H000626-1a-A39-A51 #`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
sample_normal=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sample_normal=1196553-TP6_FDHE19H000592-1a-A14-A30

IDlist=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/group/batch4_ID_group_WES_pairwise.csv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $IDlist))'
#sampleID=4824566 #`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
#sampleID=1196553

normal=${out_path}${sample_normal}/${sample_normal}.sorted.markdup.BQSR_hg38_bed.bam
tumor=${out_path}${sample_tumor}/${sample_tumor}.sorted.markdup.BQSR_hg38_bed.bam

mkdir ${work_path}/${sampleID}
outfile=${work_path}/${sampleID}

source activate reseq_env

samtools index ${tumor}
samtools index ${normal}
gatk Mutect2 -R $genome \
             -I ${tumor} -tumor $sample_tumor \
             -I ${normal} -normal $sample_normal \
             --germline-resource $ge \
             -L ${bed} \
             --af-of-alleles-not-in-resource 0.0000025  \
             --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
             -pon ${pon} \
             --f1r2-tar-gz ${outfile}/${sampleID}_f1r2.tar.gz \
             -O ${outfile}/${sampleID}_mutect2.vcf \
             --tmp-dir $tmp_dir  \
			 --native-pair-hmm-threads 8
	

gatk LearnReadOrientationModel -I ${outfile}/${sampleID}_f1r2.tar.gz \
                               -O ${outfile}/${sampleID}_read_orientation-model.tar.gz \
                               --tmp-dir $tmp_dir

gatk GetPileupSummaries -V $ge \
            -L ${bed} \
            -I ${normal} \
            -O ${outfile}/${sampleID}_n_getpileupsummaries_normal.table \
            --tmp-dir $tmp_dir

gatk GetPileupSummaries -V $ge \
            -L ${bed} \
            -I ${tumor} \
            -O ${outfile}/${sampleID}_t_getpileupsummaries_normal.table \
            --tmp-dir $tmp_dir

gatk CalculateContamination \
            -I ${outfile}/${sampleID}_t_getpileupsummaries_normal.table \
            -matched ${outfile}/${sampleID}_n_getpileupsummaries_normal.table \
            -O ${outfile}/${sampleID}_calculatecontamination.table \
            --tmp-dir $tmp_dir

gatk FilterMutectCalls -V ${outfile}/${sampleID}_mutect2.vcf \
                       -R $genome \
                       --contamination-table ${outfile}/${sampleID}_calculatecontamination.table \
                       -O ${outfile}/${sampleID}_filtered.vcf.gz \
                       --ob-priors ${outfile}/${sampleID}_read_orientation-model.tar.gz \
                       --tmp-dir $tmp_dir
source activate gatk

samtools mpileup -B -q 1 -f ${genome} ${normal} ${tumor} > ${outfile}/${sampleID}.mpileup
varscan somatic ${outfile}/${sampleID}.mpileup ${outfile}/${sampleID} --mpileup 1 -min-coverage 10 -min-var-freq 0.08 -somatic-p-value 0.05 --output-vcf 1 
#
varscan processSomatic ${outfile}/${sampleID}.snp.vcf
varscan processSomatic ${outfile}/${sampleID}.indel.vcf
#
varscan somaticFilter ${outfile}/${sampleID}.snp.Somatic.hc.vcf -indel-file ${outfile}/${sampleID}.indel.vcf -output-file ${outfile}/${sampleID}.snp.Somatic.hc.filter.vcf

bed=/data_group/cunyupeng/ynch_lungcancer/XWluad/rawData/ref/hg38_Exome.bed.gz
source activate py27
#
rm -r ${outfile}/manta
rm -r ${outfile}/strelka
mkdir ${outfile}/manta
mkdir ${outfile}/strelka
#
MANTA_INSTALL_PATH=/opt/service/miniconda3/envs/py27/share/manta-1.6.0-0
#
${MANTA_INSTALL_PATH}/bin/configManta.py \
  --normalBam ${normal} \
  --tumorBam ${tumor} \
  --referenceFasta ${genome} \
  --callRegions ${bed}  \
  --exome \
  --runDir ${outfile}/manta
#
${outfile}/manta/runWorkflow.py -j 8
#
STRELKA_INSTALL_PATH=/opt/service/miniconda3/envs/py27/share/strelka-2.9.10-0
${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
  --normalBam ${normal} \
  --tumorBam ${tumor} \
  --referenceFasta ${genome} \
  --indelCandidates ${outfile}/manta/results/variants/candidateSmallIndels.vcf.gz \
  --callRegions ${bed} \
  --exome \
  --runDir ${outfile}/strelka
#
python2 ${outfile}/strelka/runWorkflow.py -m local -j 8
