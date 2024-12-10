#!/usr/bin/bash

source activate porechop_abi_v0.5.0

t=1
s=test

ln -s ../reads.fastq.gz $s.fastq.gz

/usr/bin/time porechop_abi \
    --ab_initio \
    --verbosity 1 \
    --threads $t \
    --input $s.fastq.gz \
    --output $s.trimmed.fastq.gz \
    --temp_dir ${s}_tmp \
    >$s.porechop_abi.stdout 2>$s.porechop_abi.stderr
rm -r test_tmp
