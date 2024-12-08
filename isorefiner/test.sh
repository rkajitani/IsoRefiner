#rm -rf isorefiner_filter_work
#isorefiner filter -i ../test/bambu/out_transcript.gtf -r ../test/reads.fastq.gz -g ../test/genome_chr22.fa

rm -rf isorefiner_refine_work
isorefiner refine -i ../test/*/filter/out_transcript.gtf -r ../test/reads.fastq.gz -g ../test/genome_chr22.fa -a ../test/annot_reduced.gtf
