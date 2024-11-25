#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
    print("usage:", sys.argv[0], "in.gtf bambu_counts.tsv", file=sys.stderr)
    sys.exit(1)


def parse_gtf_attr(attr_field):
    gene_id = ""
    transcript_id = ""
    for attr_str in attr_field.split(";"):
        attr = attr_str.strip().split(" ")
        if len(attr) != 2:
            continue
        if attr[0] == "gene_id":
            gene_id = attr[1].strip('"')
        elif attr[0] == "transcript_id":
            transcript_id = attr[1].strip('"')
    return gene_id, transcript_id


positive_set = set()
with open(sys.argv[2]) as fin:
    fin.readline()
    for ln in fin:
        f = ln.rstrip("\n").split("\t")
        if float(f[2]) > 0.0:
            positive_set.add(f[0])

with open(sys.argv[1]) as fin:
    fin.readline()
    for ln in fin:
        if len(ln) == 0 or ln[1] == "#":
            print(ln, end="")
            continue

        f = ln.rstrip("\n").split("\t")
        gene_id, transcript_id = parse_gtf_attr(f[8])
        if (gene_id != "" and
            transcript_id != "" and
            transcript_id in positive_set):
            print(ln, end="")
