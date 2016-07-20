#!/bin/bash
# http://barcwiki.wi.mit.edu/wiki/SOPs/variant_calling_GATK
# The above may be useful, particularly the validate Varients step?
# This script will expect  an argument which will be fed from parallel
# like so 
out_dir=/nfs/jjmdata/spiro/
ref_genome=/nfs/zorba/jjmdata/spiro/genome.fa
gold_indels=/nfs/zorba/jjmdata/GATK/Mills_and_1000G_gold_standard.indels.b37.vcf
dbsnp=/nfs/zorba/jjmdata/GATK/dbsnp_138.b37.vcf
input=accepted_hits_karyo

for fol in `cat GATK_list.txt`; do {
cd /nfs/zorba/jjmdata/spiro/Sample_${fol}_thout/
pwd
echo "Now reordering"
java -jar $PICARD/ReorderSam.jar I=accepted_hits.bam O=${input}.bam REFERENCE=$ref_genome

} done
# The below section (through line 50) was run for each sample ID (fed as the first argument to bash ($1)
cd /nfs/zorba/jjmdata/spiro/Sample_$1_thout/
echo "Now adding read groups and sorting"
java -jar $PICARD/AddOrReplaceReadGroups.jar I=${input}.bam O=rg_added_${input}.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=$1 
echo "Now marking duplicates"
java -jar $PICARD/MarkDuplicates.jar I=rg_added_${input}.bam O=dedupped_${input}.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output_${input}.metrics 

echo "Now trimming"
java -jar $GATK -T SplitNCigarReads -R $ref_genome -I dedupped_${input}.bam -o split_${input}.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -fixNDN

echo "Now realigning about indels"
# Here run the optional steps of realignment about the indels
# http://gatkforums.broadinstitute.org/discussion/2800/howto-perform-local-realignment-around-indels
java -jar $GATK -T RealignerTargetCreator -R $ref_genome -I split_${input}.bam -known $gold_indels -U ALLOW_N_CIGAR_READS -o realignment_targets_${input}.list

java -jar $GATK -T IndelRealigner -R $ref_genome -I split_${input}.bam -targetIntervals realignment_targets_${input}.list -known $gold_indels -o realigned_reads_${input}.bam

# Here run BQSR
echo "Now running BQSR: BaseRecalibrator Part1"
# http://gatkforums.broadinstitute.org/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr
java -jar $GATK -T BaseRecalibrator -R $ref_genome -I realigned_reads_${input}.bam -knownSites $dbsnp -knownSites $gold_indels -o recal_data_${input}.table

echo "Now running BQSR: BaseRecalibrator Part 2"
java -jar $GATK -T BaseRecalibrator -R $ref_genome -I realigned_reads_${input}.bam -knownSites $dbsnp -knownSites $gold_indels -BQSR recal_data_${input}.table -o post_recal_data_${input}.table 

echo "Now running BQSR: AnalyzeCovariates"
java -jar $GATK -T AnalyzeCovariates -R $ref_genome -before recal_data_${input}.table -after post_recal_data_${input}.table -plots recalibration_plots_${input}.pdf

echo "Now running Variant Calling"
java -jar $GATK -T HaplotypeCaller -R $ref_genome -I realigned_reads_${input}.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --emitRefConfidence GVCF -o output.g.vcf

echo "Now removing extra files"
#rm dedupped_${input}.ba* ${input}.ba* rg_added_${input}.ba* split_${input}.ba*   
rm dedupped_${input}.ba* rg_added_${input}.ba* split_${input}.ba*   

# Here is the end of the per sample called script (using $1 as the first argument
# At this point each sample will have an output.g.vcf file

input=''
for i in $(cat $input_list); do {
input="$input --variant ${out_dir}Sample_${i}_thout/output.g.vcf"
} done

echo "Now running joint Variant Calling on all samples"
java -jar $GATK -T GenotypeGVCFs -R $ref_genome --dbsnp $dbsnp $input -o $out_dir/merged_output.vcf

echo "Now running Variant Filtering"
java -jar $GATK -T VariantFiltration -R $ref_genome -V merged_output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o merged_output_filtered.vcf

# Here run snpEff (to annotated the merged file resulting merged_output_filtered.ann.vcf) and snpSift (dbsnp and case control resulting in merged_output_filtered_annotated_SUI.ann.vcf and merged_output_filtered_annotated_MDD.ann.vcf) as explained in
# http://snpeff.sourceforge.net  

for a in SUI MDD; do {
file=merged_output_filtered_annotated_$a.ann.vcf
for thr in 0.05 0.0001; do { # the 0.05 will be for top34 genes, the 0.0001 will be for whole-exome
# here first filter to only include variants with subthreshold differences in case vs. control
cat $file | java -jar $snpSift filter "(CC_DOM_$a < $thr | CC_REC_$a < $thr | CC_ALL_$a < $thr | CC_GENO_$a < $thr | CC_TREND_$a < $thr)" > subthresh_thr${thr}_$file 

# Here extract relevant fields from filtered file 
java -jar $snpSift extractFields -s "," -e "." subthresh_thr${thr}_$file CHROM POS RS REF ALT Cases_$a Controls_$a CC_DOM_$a CC_REC_$a CC_ALL_$a CC_GENO_$a CC_TREND_$a "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].EFFECT" > extracted_thr${thr}_fields_$a.txt
} done

# here to select the top34 genes from the file thresholded at p<0.05 uncorrected
thr=0.05
head extracted_thr${thr}_fields_$a.txt -n 1 > subthresh_pthr${thr}_genes_extracted_fields_$a.txt
for i in $(cat top34_genes.txt); do {
cat subthresh_extracted_fields_SUI.txt | grep $i >> top34_subthresh_genes_extracted_fields_SUI.txt
cat subthresh_extracted_fields_MDD.txt | grep $i >> top34_subthresh_genes_extracted_fields_MDD.txt
} done
} done



