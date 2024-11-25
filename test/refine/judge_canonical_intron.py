#!/usr/bin/env python

import re
import sys


if len(sys.argv) != 2:
    print("usage:", sys.argv[0], "intron_edge.tsv", file=sys.stderr)
    sys.exit(1)

#canonical_set = set([
#    ("GT", "AG"),
#    ("GC", "AG"),
#    ("AT", "AC")
#])
canonical_set = set([("GT", "AG")])

exon_id_re = re.compile(r'^(\S+),\d+$')
id_set = set()
non_canonical_id_set = set()
with open(sys.argv[1]) as fin:
    for ln in fin:
        f = ln.rstrip("\n").split("\t")
        m = exon_id_re.match(f[0])
        if not m:
            continue
        transcript_id = m.group(1)

        id_set.add(transcript_id)
        if (f[1], f[2]) not in canonical_set:
            non_canonical_id_set.add(transcript_id)

print("transcript_id", "non_canonical_flag", sep="\t")
for transcript_id in id_set:
    non_canonical_flag = 1
    if transcript_id not in non_canonical_id_set:
        non_canonical_flag = 0
    print(transcript_id, non_canonical_flag, sep="\t")
