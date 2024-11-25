#!/usr/bin/bash

source activate espresso_v1.5.0

t=4
start_dir=$PWD
bam_dir=../map_genome

samples=(test)

rm samples.tsv &>/dev/null
for s in ${samples[@]}
do
    echo -e "$bam_dir/$s.sorted.bam\t$s" >>samples.tsv
done

ln -s ../genome_chr22.fa genome.fa
ln -s ../annot_reduced.gtf annot.gtf

/usr/bin/time perl $(which ESPRESSO_S.pl) -T $t -L samples.tsv -F genome.fa -A annot.gtf -O espresso_work >ESPRESSO_S.stdout 2>ESPRESSO_S.stderr
for i in $(seq 0 $(expr ${#samples[@]} - 1))
do
    /usr/bin/time perl $(which ESPRESSO_C.pl) -X $i -T $t -I espresso_work -F genome.fa >ESPRESSO_C_$i.stdout 2>ESPRESSO_C_$i.stderr
done
/usr/bin/time perl $(which ESPRESSO_Q.pl) -T $t -L espresso_work/samples.tsv.updated -A annot.gtf >ESPRESSO_Q.stdout 2>ESPRESSO_Q.stderr
ln -s espresso_work/samples_N2_R0_updated.gtf out_transcript.gtf
