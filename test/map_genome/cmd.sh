#!/usr/bin/bash

t=4
s=test

ln -s ../genome_chr22.fa genome.fa
ln -s ../porechop_abi/test.trimmed.fastq.gz $s.fastq.gz

/usr/bin/time minimap2 \
    -t $t \
    -a \
    -x splice \
    -ub \
    -k14 \
    --secondary=no \
    genome.fa \
    $s.fastq.gz \
    2> $s.stderr \
    | samtools view -Sb > $s.bam
 
/usr/bin/time samtools sort -O BAM -@ $t -m 2G $s.bam >$s.sorted.bam 2>$s.sort.stderr
samtools index $s.sorted.bam
