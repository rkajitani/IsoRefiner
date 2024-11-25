#!/usr/bin/env python

import re
import sys
from collections import defaultdict


if len(sys.argv) != 2:
    print("usage:", sys.argv[0], "annot.gtf", file=sys.stderr)
    sys.exit(1)

transcript_re = re.compile(r'transcript_id "([^"]*)"')
n_exon_dict = defaultdict(int)
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
        n_exon_dict[transcript_id] += 1
        print(f[0], int(f[3]) - 1, f[4], f"{transcript_id},{n_exon_dict[transcript_id]}", ".", f[6], sep="\t")
