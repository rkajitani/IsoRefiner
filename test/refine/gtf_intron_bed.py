#!/usr/bin/env python

import re
import sys
from collections import defaultdict


if len(sys.argv) != 2:
    print("usage:", sys.argv[0], "annot.gtf", file=sys.stderr)
    sys.exit(1)


class ExonPos:
    def __init__(self):
        self.seq_name = ""
        self.strand = ""
        self.pos_list = []

    def add(self, seq_name, start, end, strand):
        self.seq_name = seq_name
        self.pos_list.append(start)
        self.pos_list.append(end)
        self.strand = strand

    def get_splice_sites(self):
        if len(self.pos_list) < 4:
            return []
        else:
            return sorted(self.pos_list)[1:-1]

        
transcript_re = re.compile(r'transcript_id "([^"]*)"')
exon_pos_dict = defaultdict(ExonPos)
with open(sys.argv[1]) as fin:
    for ln in fin:
        if len(ln) == 0 or ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if f[2] != "exon":
            continue

        m = transcript_re.search(f[8])
        if not m:
            continue
        transcript_id = m.group(1)
        exon_pos_dict[transcript_id].add(f[0], int(f[3]) - 1, int(f[4]), f[6])

for transcript_id, exon_pos in exon_pos_dict.items():
    splice_site_list = exon_pos.get_splice_sites()
    if len(splice_site_list) < 2:
        continue
    n_exon = 0
    for i in range(0, len(splice_site_list), 2):
        n_exon += 1
        print(
            exon_pos.seq_name,
            splice_site_list[i],
            splice_site_list[i + 1],
            f"{transcript_id},{n_exon}",
            ".",
            exon_pos.strand,
            sep="\t"
        )
