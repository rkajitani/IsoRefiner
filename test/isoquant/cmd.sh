#!/usr/bin/bash

source activate isoquant_v3.3.1

t=4
bam_dir=../map_genome

samples=(test)

echo "#isoquant" >input_bam.txt
for s in ${samples[@]}
do
    echo "$bam_dir/$s.sorted.bam:$s" >>input_bam.txt
done

ln -s ../genome_chr22.fa genome.fa
ln -s ../annot_reduced.gtf annot.gtf


/usr/bin/time isoquant.py \
    --threads $t \
    --reference genome.fa \
    --genedb annot.gtf \
    --complete_genedb \
    --bam_list input_bam.txt \
    --data_type nanopore \
    --stranded none \
    --output isoquant_out \
    --transcript_quantification unique_only \
    --gene_quantification unique_only \
    --matching_strategy default \
    --splice_correction_strategy default_ont \
    --model_construction_strategy default_ont \
    --no_secondary \
    --check_canonical \
    --count_exons \
    >isoquant.stdout 2>isoquant.stderr

ln -s isoquant_out/isoquant/isoquant.transcript_models.gtf out_transcript.gtf
