#!/usr/bin/env python

import sys
import polars as pl


if len(sys.argv) != 2:
    print("usage:", sys.argv[0], "transcript_info.tsv > result.tsv", file=sys.stderr)
    sys.exit(1)

dist_th = 20

df = pl.read_csv(sys.argv[1], separator="\t")
print("\t".join(df.columns))
for row in df.rows(named=True):
    transcript_id = row["transcript_id"]
    counter_id = row["counter_id"]
    intron_pos_dist = row["intron_pos_dist"]
    n_intron = row["n_intron"]
    counter_n_intron = row["counter_n_intron"]
    transcript_depth = row["transcript_depth"]
    counter_depth = row["counter_depth"]
    non_canonical_flag = row["non_canonical_flag"]
    counter_non_canonical_flag = row["counter_non_canonical_flag"]
    transcript_len = row["transcript_len"]
    counter_len = row["counter_len"]

    target_flag = False
    
    if intron_pos_dist == 0:
        if (n_intron <= counter_n_intron and transcript_len <= counter_len and transcript_depth <= counter_depth and
            (n_intron < counter_n_intron or transcript_len < counter_len or transcript_depth < counter_depth)):
                target_flag = True
    elif (intron_pos_dist <= dist_th and 
        non_canonical_flag == 1 and
        counter_non_canonical_flag == 0 and
        n_intron <= counter_n_intron and
        transcript_len <= counter_len and
        transcript_depth <= counter_depth):
        target_flag = True

    if target_flag:
        print("\t".join([str(_) for _ in row.values()]))
