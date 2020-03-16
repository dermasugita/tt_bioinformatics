#!/bin/bash
set -eu
NR1=$1
NR2=$2
TR1=$3
TR2=$4
reference='../reference/chr11.fa'

touch ../tmp/log

bwa mem -M -t 2 $reference $NR1 $NR2 | samtools view -bS - > ../tmp/N-WEX.bam
bwa mem -M -t 2 $reference $TR1 $TR2 | samtools view -bS - > ../tmp/T-WEX.bam

echo 'bwa mem finished' >> ../tmp/log

gatk AddOrReplaceReadGroups -I ../tmp/N-WEX.bam -O ../tmp/N-WEX.renamed.bam -LB N-WEX -PL illumina -PU N-WEX -SM N-WEX  
gatk AddOrReplaceReadGroups -I ../tmp/T-WEX.bam -O ../tmp/T-WEX.renamed.bam -LB T-WEX -PL illumina -PU T-WEX -SM T-WEX  

echo 'AddOrReplaceReadGroups finished' >> ../tmp/log 
gatk MarkDuplicatesSpark -I ../tmp/N-WEX.renamed.bam -O ../tmp/N-WEX.marked.bam -M ../tmp/N-WEX.marked.txt  
gatk MarkDuplicatesSpark -I ../tmp/T-WEX.renamed.bam -O ../tmp/T-WEX.marked.bam -M ../tmp/T-WEX.marked.txt  

echo 'AddOrReplaceReadGroups finished' >> ../tmp/log 
echo 'MarkDupicatesSpark finished' >> ../tmp/log

gatk BaseRecalibrator -I ../tmp/N-WEX.marked.bam -R $reference --known-sites ../dbSNP/dbsnp_146.hg38.chr11.vcf.gz --known-sites ../dbSNP/Mills_and_1000G_gold_standard.indels.hg38.chr11.vcf.gz -O ../tmp/N_bqsr.table 
gatk BaseRecalibrator -I ../tmp/T-WEX.marked.bam -R $reference --known-sites ../dbSNP/dbsnp_146.hg38.chr11.vcf.gz --known-sites ../dbSNP/Mills_and_1000G_gold_standard.indels.hg38.chr11.vcf.gz -O ../tmp/T_bqsr.table 

echo 'BaseRecalibrator finished' >> ../tmp/log

gatk ApplyBQSR -R $reference -I ../tmp/N-WEX.marked.bam --bqsr-recal-file ../tmp/N_bqsr.table -O ../bam/N-WEX.marked.bqsr.bam 
gatk ApplyBQSR -R $reference -I ../tmp/T-WEX.marked.bam --bqsr-recal-file ../tmp/T_bqsr.table -O ../bam/T-WEX.marked.bqsr.bam 

echo 'ApplyBQSR finished' >> ../tmp/log

gatk Mutect2 -R $reference -I ../bam/N-WEX.marked.bqsr.bam -I ../bam/T-WEX.marked.bqsr.bam \
	-normal N-WEX -tumor T-WEX --germline-resource ../dbSNP/af-only-gnomad.hg38.chr11.vcf.gz \
	--f1r2-tar-gz ../tmp/f1r2.tar.gz -O ../tmp/unfiltered.vcf 

echo 'Finished Mutect2' >> ../tmp/log

gatk LearnReadOrientationModel -I ../tmp/f1r2.tar.gz -O ../tmp/read-orientation-model.tar.gz 

echo 'Finished LearnReadOrientationModel' >> ../tmp/log

gatk GetPileupSummaries -I ../bam/N-WEX.marked.bqsr.bam -V ../dbSNP/small_exac_common_3.hg38.chr11.vcf.gz -L ../dbSNP/small_exac_common_3.hg38.chr11.vcf.gz -O ../tmp/N_pileup.table 
gatk GetPileupSummaries -I ../bam/T-WEX.marked.bqsr.bam -V ../dbSNP/small_exac_common_3.hg38.chr11.vcf.gz -L ../dbSNP/small_exac_common_3.hg38.chr11.vcf.gz -O ../tmp/T_pileup.table 

echo 'Finished GetPileupSummaries' >> ../tmp/log

gatk CalculateContamination -I ../tmp/T_pileup.table -matched ../tmp/N_pileup.table -tumor-segmentation ../tmp/segments.table -O ../tmp/contamination.table 

echo 'Finished CalculateContamination' >> ../tmp/log

gatk FilterMutectCalls -V ../tmp/unfiltered.vcf -R $reference --tumor-segmentation ../tmp/segments.table --contamination-table ../tmp/contamination.table \
	--ob-priors ../tmp/read-orientation-model.tar.gz -O ../vcf/filtered.vcf
