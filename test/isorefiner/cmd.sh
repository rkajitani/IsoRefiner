#!/usr/bin/bash

reads=../reads.fastq.gz
genome=../genome_chr22.fa
ref_annot=../annot_reduced.gtf
n_threads=1

# Clean results before test 
rm -r isorefiner_refined.gtf isorefiner_trans_struct_wf_work  &>/dev/null

# Run the workflow.
# Final result: isorefiner_refined.gtf
isorefiner trans_struct_wf -r $reads -g $genome -a $ref_annot -t $n_threads
