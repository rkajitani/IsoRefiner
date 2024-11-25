#!/usr/bin/bash

t=16
sqanti_sim=/mnt/fsx/user/kajitt3k/tools/SQANTI-SIM_v0.2.0/sqanti-sim.py
start_dir=$PWD

./gtf_filt_zero_count.py bambu_out/extended_annotations.gtf bambu_out/counts_transcript.txt >out_transcript.gtf


for m in 3 10 20 30
do
	mkdir -p min_support_$m
	cd min_support_$m

	cat <<_EOS >cmd.sh
#!/usr/bin/bash

source activate sqanti_sim_v0.2.0

t=$t
sqanti_sim=$sqanti_sim
m=$m

_EOS

	cat <<'_EOS' >>cmd.sh
/usr/bin/time python $sqanti_sim eval \
    --transcriptome ../out_transcript.gtf \
    --gtf ../annot.gtf \
    --genome ../genome.fa \
    -i ../../sim/sqanti-sim_index.tsv \
    -o sqanti_sim_eval \
	--min_support $m \
    -k $t \
    >sqanti_sim_eval.stdout 2>sqanti_sim_eval.stderr
_EOS

	bash cmd.sh &

	cd $start_dir
done
wait
