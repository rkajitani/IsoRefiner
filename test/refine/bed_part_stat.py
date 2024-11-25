#!/usr/bin/env python

import re
import sys
from collections import defaultdict
from operator import itemgetter


class Transcript:
    def __init__(self):
        self.n = 0
        self.length = 0

    def add(self, length):
        self.n += 1
        self.length += length



if len(sys.argv) != 2:
    print("usage:", sys.argv[0], "parts.bed", file=sys.stderr)
    sys.exit(1)

part_id_re = re.compile(r'^(\S+),\d+$')
part_stat_dict = defaultdict(Transcript)
with open(sys.argv[1]) as fin:
    for ln in fin:
        f = ln.rstrip("\n").split("\t")

        m = part_id_re.match(f[3])
        if not m:
            continue
        transcript_id = m.group(1)
        part_stat_dict[transcript_id].add(int(f[2]) - int(f[1]))

print("transcript_id", "n_parts", "sum_len", sep="\t")
for transcript_id, part_stat in part_stat_dict.items():
    print(transcript_id, part_stat.n, part_stat.length, sep="\t")
