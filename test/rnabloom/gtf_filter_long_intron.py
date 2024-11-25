#!/usr/bin/env python

import re
import sys
from collections import defaultdict


class ExonPos:
    def __init__(self):
        self.seq_name = ""
        self.pos_list = []

    def add(self, seq_name, start, end):
        self.seq_name = seq_name
        self.pos_list.append(start)
        self.pos_list.append(end)

    def get_splice_sites(self):
        if len(self.pos_list) < 4:
            return []
        else:
            return sorted(self.pos_list)[1:-1]


if len(sys.argv) != 3:
    print("usage:", sys.argv[0], "annot.gtf max_intron_len", file=sys.stderr)
    sys.exit(1)


max_intron_len = int(sys.argv[2])
        
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
        exon_pos_dict[transcript_id].add(f[0], int(f[3]) - 1, int(f[4]))

rm_set = set()
for transcript_id, exon_pos in exon_pos_dict.items():
    splice_site_list = exon_pos.get_splice_sites()
    if len(splice_site_list) < 2:
        continue
    for i in range(0, len(splice_site_list), 2):
        if splice_site_list[i + 1] - splice_site_list[i] > max_intron_len:
            rm_set.add(transcript_id)

with open(sys.argv[1]) as fin:
    for ln in fin:
        if len(ln) == 0 or ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")

        m = transcript_re.search(f[8])
        if not m:
            continue
        transcript_id = m.group(1)

        if transcript_id not in rm_set:
            print(ln, end="")
