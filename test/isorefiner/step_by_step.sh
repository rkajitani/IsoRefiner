#!/usr/bin/bash

reads=../reads.fastq.gz
genome=../genome_chr22.fa
ref_annot=../annot_reduced.gtf


# Clean results before test 
rm -r isorefiner_* &>/dev/null


isorefiner trim -r $reads
isorefiner map -r isorefiner_trimmed.fastq.gz -g $genome

# Run Mapping-based tools. The loop below can be parallelized.
for tool in bambu espresso isoquant stringtie
do
	isorefiner run_$tool -b isorefiner_mapped.bam -g $genome -a $ref_annot
	isorefiner filter -i isorefiner_$tool.gtf -r isorefiner_trimmed.fastq.gz -g $genome -o isorefiner_filtered_$tool.gtf -d isorefiner_filter_${tool}_work
done

# Run de novo assembly tool
tool=rnabloom
isorefiner run_$tool -r isorefiner_trimmed.fastq.gz -g $genome
isorefiner filter -i isorefiner_$tool.gtf -r isorefiner_trimmed.fastq.gz -g $genome -o isorefiner_filtered_$tool.gtf -d isorefiner_filter_${tool}_work

isorefiner refine -i isorefiner_filtered_*.gtf -r isorefiner_trimmed.fastq.gz -g $genome -a $ref_annot

# Final result: isorefiner_refined.gtf
