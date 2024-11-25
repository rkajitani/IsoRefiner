#!/usr/bin/bash

source activate bambu_v3.4.0

t=4
bam_dir=../map_genome
sample_list=(test)

for s in ${sample_list[@]}
do
    bam_list+=($bam_dir/$s.sorted.bam)
done

ln -s ../genome_chr22.fa genome.fa
ln -s ../annot_reduced.gtf annot.gtf

/usr/bin/time Rscript --slave --vanilla bambu.R $t genome.fa annot.gtf ${bam_list[@]} >bambu.stdout 2>bambu.stderr
./gtf_filt_zero_count.py bambu_out/extended_annotations.gtf bambu_out/counts_transcript.txt >out_transcript.gtf
