#rm -r isorefiner_trim_work
#isorefiner trim -r ../test/reads.fastq.gz
#
#rm -rf isorefiner_filter_work
#isorefiner filter -i ../test/stringtie/out.gtf -r ../test/reads.fastq.gz -g ../test/genome_chr22.fa
#
#rm -r isorefiner_map_work
#isorefiner map -r ../test/reads.fastq.gz -g ../test/genome_chr22.fa
#
#rm -rf isorefiner_bambu_work
#isorefiner run_bambu -b  ../test/map_genome/test.sorted.bam -g ../test/genome_chr22.fa -a ../test/annot_reduced.gtf
#
#rm -rf isorefiner_espresso_work
#isorefiner run_espresso -b  ../test/map_genome/test.sorted.bam -g ../test/genome_chr22.fa -a ../test/annot_reduced.gtf
#
#rm -rf isorefiner_isoquant_work
#isorefiner run_isoquant -b  ../test/map_genome/test.sorted.bam -g ../test/genome_chr22.fa -a ../test/annot_reduced.gtf
#
#rm -rf isorefiner_stringtie_work
#isorefiner run_stringtie -b  ../test/map_genome/test.sorted.bam -g ../test/genome_chr22.fa -a ../test/annot_reduced.gtf

rm -rf isorefiner_rnabloom_work
isorefiner run_rnabloom -r ../test/reads.fastq.gz -g ../test/genome_chr22.fa --max_mem 4g

#rm -rf isorefiner_refine_work
#isorefiner refine -i ../test/*/filter/out_transcript.gtf -r ../test/reads.fastq.gz -g ../test/genome_chr22.fa -a ../test/annot_reduced.gtf
