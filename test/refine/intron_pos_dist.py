#!/usr/bin/env python

import re
import sys
from collections import defaultdict
from operator import itemgetter


class Overlap:
    def __init__(self):
        self.n = 0
        self.dist = 0
    def calc_dist(self, start1, end1, start2, end2):
        self.dist += abs(start1 - start2) + abs(end1 - end2)
        self.n += 1


if len(sys.argv) != 3:
    print("usage:", sys.argv[0], "bedtools_intersect.gtf intron_n_len.tsv", file=sys.stderr)
    sys.exit(1)

n_intron_dict = dict()
with open(sys.argv[2]) as fin:
    fin.readline()
    for ln in fin:
        f = ln.rstrip("\n").split("\t")
        n_intron_dict[f[0]] = int(f[1])

exon_id_re = re.compile(r'^(\S+),\d+$')
ovl_dict = defaultdict(lambda: defaultdict(Overlap))
with open(sys.argv[1]) as fin:
    for ln in fin:
        f = ln.rstrip("\n").split("\t")

        m = exon_id_re.match(f[3])
        if not m:
            continue
        transcript_id = m.group(1)

        m = exon_id_re.match(f[9])
        if not m:
            continue
        counter_id = m.group(1)

        if n_intron_dict[transcript_id] <= n_intron_dict[counter_id]:
            ovl_dict[transcript_id][counter_id].calc_dist(
                int(f[1]),
                int(f[2]),
                int(f[7]),
                int(f[8])
            )

print("transcript_id", "counter_id", "intron_pos_dist", "n_intron", "counter_n_intron", sep="\t")
for transcript_id, each_ovl_dict in ovl_dict.items():
    n_intron = n_intron_dict[transcript_id]
    min_dist_id = ""
    min_dist = 0
    for counter_id, ovl in each_ovl_dict.items():
        if n_intron == ovl.n and (min_dist_id == "" or min_dist > ovl.dist):
            min_dist = ovl.dist
            min_dist_id = counter_id
    if min_dist_id != "":
        print(transcript_id, min_dist_id, min_dist, n_intron, n_intron_dict[counter_id], sep="\t")
