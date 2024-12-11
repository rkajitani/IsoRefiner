#!/usr/bin/bash

t=4
max_mem=400g
export JAVA_TOOL_OPTIONS="-Xmx$max_mem"
sample_list=(test)

ln -s ../porechop_abi/test.trimmed.fastq.gz test.fastq.gz

for s in ${sample_list[@]}
do
    fastq_list+=($s.fastq.gz)
done

/usr/bin/time rnabloom -t $t -outdir rnabloom_out -long ${fastq_list[@]} >rnabloom.stdout 2>rnabloom.stderr


t=4
db_dir=.
db_name=genome_gmap
min_cov=0.5
min_idt=0.95
max_intron_len=100000
genome_fa=../genome_chr22.fa
transcript_fa=./rnabloom_out/rnabloom.transcripts.fa 

/usr/bin/time gmap_build -D $db_dir -d $db_name $genome_fa >gmap_build.stdout 2>>gmap_build.stderr

/usr/bin/time gmap \
    -D $db_dir \
    -d $db_name \
    -t $t \
    -f 2 \
    -n 1 \
    --no-chimeras \
    --min-trimmed-coverage=$min_cov \
    --min-identity=$min_idt \
    --max-intronlength-middle=$max_intron_len \
    $transcript_fa >gmap.gff3 2>gmap.stderr

gffread -E -T -o gmap.gtf gmap.gff3
./gtf_filter_long_intron.py gmap.gtf $max_intron_len >gmap_filt.gtf
rm -r genome_gmap
ln -s gmap_filt.gtf out_transcript.gtf
