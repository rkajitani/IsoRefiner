#!/usr/bin/env python

import sys
import pysam


if len(sys.argv) != 6:
    print("usage:", sys.argv[0], "input.bam output_bam max_indel max_clip min_idt(0-1)")
    sys.exit(1)

in_bam = sys.argv[1]
out_bam = sys.argv[2]
max_indel = int(sys.argv[3])
max_clip = int(sys.argv[4])
min_idt = float(sys.argv[5])

with pysam.AlignmentFile(in_bam, "rb") as fin, \
     pysam.AlignmentFile(out_bam, "wb", header=fin.header) as fout:
        for aln in fin:
            valid_flag = True

            idt = 1.0 - (aln.get_tag("NM") / aln.query_alignment_length)
            if idt < min_idt:
                continue

            for cig_ope, cig_len in aln.cigartuples:
                if (((cig_ope == 1 or cig_ope == 2) and cig_len > max_indel) or
                    (cig_ope == 4 or cig_ope == 5) and cig_len > max_clip):
                    valid_flag = False
                    break

            if valid_flag:
                fout.write(aln)
