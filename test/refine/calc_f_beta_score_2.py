#!/usr/bin/env python

import sys
import polars as pl

df = pl.read_csv(sys.argv[1], separator="\t")

filt_list = [
    ("transcript", "count10_removed", "recall"),
    ("transcript", "all_ref", "precision"),
]

filt_df_list = list()
for filt in filt_list:
    filt_df_list.append((df
        .filter(
            (pl.col("level") == filt[0]) & 
            (pl.col("condition") == filt[1])
        )
        .select(["tool", filt[2]])
    ))


def f_beta(struct):
    beta = 1
    r = struct["recall"]
    p = struct["precision"]
    return (beta**2 + 1) * (p * r) / (((beta**2) * p) + r)

joined_df = (filt_df_list[0]
    .join(filt_df_list[1], on="tool", how="inner")
    .with_columns(pl.struct(["recall", "precision"]).map_elements(f_beta).alias("f_beta_score"))
    .sort(by="f_beta_score", descending=True)
    .with_row_count("f_beta_score_rank", offset=1)
    .select(["tool", "recall", "precision", "f_beta_score", "f_beta_score_rank"])
)


joined_df.write_csv(sys.stdout, separator="\t")
