#!/usr/bin/bash

source activate stringtie_v2.2.1

t=4
bam_dir=../map_genome
sample_list=(test)

for s in ${sample_list[@]}
do
    bam_list+=($bam_dir/$s.sorted.bam)
done

ln -s ../annot_reduced.gtf annot.gtf

/usr/bin/time stringtie -o out.gtf -G annot.gtf -L ${bam_list[@]} >stringtie.stdout 2>stringtie.stderr
perl -F"\t" -ane '$F[6] = "+" if (not(/^#/) and $F[6] eq "."); print(join("\t", @F))' out.gtf >out_transcript.gtf
