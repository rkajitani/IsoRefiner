#!/usr/bin/bash

### Merge
ln -s ../annot_reduced.gtf ref.gtf

raw_gtf_files=(\
    ../bambu/filter/out_transcript.gtf \
    ../espresso/filter/out_transcript.gtf \
    ../isoquant/filter/out_transcript.gtf \
    ../rnabloom/filter/out_transcript.gtf \
    ../stringtie/filter/out_transcript.gtf \
)

for i in $(seq 0 $(expr ${#raw_gtf_files[@]} - 1))
do
    ln -s ${raw_gtf_files[$i]} set_$i.gtf
    gtf_files+=(set_$i.gtf)
done

/usr/bin/time gffcompare -o merge -p cons -r ref.gtf ref.gtf ${gtf_files[@]} >gffcompare_merge.stdout 2>gffcompare_merge.stdout


### Strand correction
/usr/bin/time gffcompare -r ref.gtf merge.combined.gtf -o flip >gffcompare_flip.log 2>gffcompare_flip.stderr
sed 1d flip.merge.combined.gtf.tmap | perl -ane 'print($F[4], "\n") if ($F[2] eq "s")' >opposite_strand.list
./gtf_flip_strand.py merge.combined.gtf opposite_strand.list >flipped.gtf 2>gtf_flip_strand.stderr


### Correction based on intron-distance
t=4
max_indel=20
max_clip=200
min_idt=0.9
min_cov=95
min_mean_depth=1

ln -s ../flip/flipped.gtf raw.gtf
ln -s ../genome_chr22.fa genome.fa
ln -s ../reads.fastq.gz reads.fastq.gz

gffread -w asm.fa -g genome.fa flipped.gtf
/usr/bin/time minimap2 -ax map-ont --secondary=no -t $t asm.fa reads.fastq.gz 2>mm2.log | samtools view -b -F 2308 - >raw.bam
samtools sort -@ $t -m 2G raw.bam >sorted.bam 2>sort.log
./bam_filter.py sorted.bam filt.bam $max_indel $max_clip $min_idt
samtools index filt.bam
samtools coverage filt.bam >cov_raw.tsv
df_select.py cov_raw.tsv '#rname' coverage meandepth | perl -pne 's/#rname/qry_id/' >cov.tsv

./gtf_exon_bed.py flipped.gtf >exon.bed
./bed_part_stat.py exon.bed >exon_n_len.tsv

./gtf_intron_bed.py flipped.gtf >intron.bed
./bed_part_stat.py intron.bed >intron_n_len.tsv
(bedtools intersect -s -wao -a intron.bed -b intron.bed | perl -ane 'print if ($F[3] ne $F[9])') >intron_self_ovl.tsv 2>bedtools_intersect.log
./intron_pos_dist.py intron_self_ovl.tsv intron_n_len.tsv >intron_pos_dist.tsv
bedtools getfasta -nameOnly -s -bed intron.bed -fi genome.fa \
    | perl -pne 's/\([+-]\)//' \
    | seqkit fx2tab \
    | perl -ane 'print(join("\t", ($F[0], substr($F[1], 0, 2), substr($F[1], -2, 2))), "\n")' \
    >intron_edge.tsv
./judge_canonical_intron.py intron_edge.tsv >non_canonical_flag.tsv

tsv_left_join.py -d -l 1 intron_pos_dist.tsv <(cut -f1,3 cov.tsv | perl -pne 's/qry_id\tmeandepth/counter_id\tcounter_depth/ if ($. == 1)') >tmp
tsv_left_join.py -d -l 1 tmp <(cut -f1,3 cov.tsv | perl -pne 's/qry_id\tmeandepth/transcript_id\ttranscript_depth/ if ($. == 1)') >tmp2
tsv_left_join.py -d -l 1 tmp2 <(perl -pne 's/non_canonical_flag/counter_non_canonical_flag/ if ($. == 1)' non_canonical_flag.tsv) >tmp3
tsv_left_join.py -d -l 1 tmp3 non_canonical_flag.tsv >tmp4
tsv_left_join.py -d -l 1 tmp4 <(cut -f1,3 exon_n_len.tsv | perl -pne 's/sum_len/counter_len/ if ($. == 1)') >tmp5
tsv_left_join.py -d -l 1 tmp5 <(cut -f1,3 exon_n_len.tsv | perl -pne 's/sum_len/transcript_len/ if ($. == 1)') >tmp6
df_select.py tmp6 \
    transcript_id \
    counter_id \
    intron_pos_dist \
    n_intron \
    counter_n_intron \
    transcript_depth \
    counter_depth \
    non_canonical_flag \
    counter_non_canonical_flag \
    transcript_len \
    counter_len \
    >intron_joined.tsv
rm tmp tmp2 tmp3 tmp4 tmp5 tmp6

./filt_by_intron_info.py intron_joined.tsv >artefact_cand.tsv
sed 1d artefact_cand.tsv | cut -f1 >artefact_cand.list
./gtf_transcript_id_grep_inv.py artefact_cand.list flipped.gtf >intron_dist_filt.gtf


### Reformat
gffcompare --strict-match -e 0 -r ref.gtf intron_dist_filt.gtf -o reform >gffcompare_reform.stdout 2>gffcompare_reform.stderr
./reform_cons_gtf.py intron_dist_filt.gtf reform.intron_dist_filt.gtf.tmap ref.gtf >reform.gtf 
