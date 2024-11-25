#!/usr/bin/env python

import re
import sys


if len(sys.argv) != 3:
    print("usage:", sys.argv[0], "transcript_id_list.txt annot.gtf", file=sys.stderr)
    sys.exit(1)

target_set = set()
with open(sys.argv[1]) as fin:
    for ln in fin:
        target_set.add(ln.rstrip("\n"))

transcript_re = re.compile(r'transcript_id "([^"]*)"')
with open(sys.argv[2]) as fin:
    for ln in fin:
        if len(ln) == 0 or ln[0] == "#":
            continue

        f = ln.rstrip("\n").split("\t")
        m = transcript_re.search(f[8])
        if (not m) or (m.group(1) not in target_set):
            print(ln, end="")
