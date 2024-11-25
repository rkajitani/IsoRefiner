#!/usr/bin/bash

t=4
max_indel=20
max_clip=200
min_idt=0.9
min_cov=95
min_mean_depth=1

ln -s ../out_transcript.gtf raw.gtf
ln -s ../../reads.fastq.gz
ln -s ../../genome_chr22.fa genome.fa

gffread -w asm.fa -g genome.fa raw.gtf
/usr/bin/time minimap2 -ax map-ont --secondary=no -t $t asm.fa reads.fastq.gz 2>mm2.log | samtools view -b -F 2308 - >raw.bam
samtools sort -@ $t -m 2G raw.bam >sorted.bam 2>sort.log
./bam_filter.py sorted.bam filt.bam $max_indel $max_clip $min_idt
samtools index filt.bam
samtools coverage filt.bam >cov.tsv

df_select.py cov.tsv "#rname" coverage meandepth | sed 1d | perl -ane "print(\$F[0], \"\n\") if (\$F[1] >= $min_cov and \$F[2] >= $min_mean_depth)" >cov_depth_filt.list
./gtf_extract_transcript.py raw.gtf cov_depth_filt.list >cov_depth_filt.gtf
ln -s cov_depth_filt.gtf out_transcript.gtf
