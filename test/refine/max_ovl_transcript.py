#!/usr/bin/env python

import re
import sys
from collections import defaultdict
from operator import itemgetter


if len(sys.argv) != 2:
    print("usage:", sys.argv[0], "bedtools_intersect.tsv", file=sys.stderr)
    sys.exit(1)

exon_id_re = re.compile(r'^(\S+),\d+$')
ovl_len_dict = defaultdict(lambda: defaultdict(int))
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

        ovl_len_dict[transcript_id][counter_id] += int(f[12])

print("transcript_id", "counter_id", "ovl_len", sep="\t")
for transcript_id, each_ovl_len_dict in ovl_len_dict.items():
    max_counter_id, max_ovl_len = max(each_ovl_len_dict.items(), key=itemgetter(1))
    print(transcript_id, max_counter_id, max_ovl_len, sep="\t")
