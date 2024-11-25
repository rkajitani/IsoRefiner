#!/usr/bin/env python

import re
import sys


if len(sys.argv) != 3:
    print("usage:", sys.argv[0], "annot.gtf transcript.list", file=sys.stderr)
    sys.exit(1)

with open(sys.argv[2]) as fin:
    transcript_set = set(fin.read().splitlines())

gene_set = set()
gene_re = re.compile(r'gene_id "([^"]*)"')
transcript_re = re.compile(r'transcript_id "([^"]*)"')
with open(sys.argv[1]) as fin:
    for ln in fin:
        if len(ln) == 0 or ln[0] == "#":
            continue

        f = ln.rstrip("\n").split("\t")
        if f[2] != "transcript":
            continue

        m = gene_re.search(f[8])
        if not m:
            continue
        gene_id = m.group(1)

        m = transcript_re.search(f[8])
        if not m:
            continue
        transcript_id = m.group(1)

        if transcript_id in transcript_set:
            gene_set.add(gene_id)

with open(sys.argv[1]) as fin:
    for ln in fin:
        if len(ln) == 0 or ln[0] == "#":
            print(ln, end="")
            continue

        f = ln.rstrip("\n").split("\t")
        if f[2] == "gene":
            m = gene_re.search(f[8])
            if not m:
                continue
            gene_id = m.group(1)
            if gene_id in gene_set:
                print(ln, end="")
        else:
            m = transcript_re.search(f[8])
            if not m:
                continue
            transcript_id = m.group(1)
            if transcript_id in transcript_set:
                print(ln, end="")
