#!/usr/bin/env python

import csv
import re
import sys
from collections import defaultdict


class GeneInfo:
    def __init__(self, seq_id, start, end, strand):
        self.seq_id = seq_id
        self.start = start
        self.end = end
        self.strand = strand 


def parse_attr(attr_field):
    attr_dict = dict()
    for attr_str in attr_field.split(";"):
        m = re.match(r'(\S+)\s+"(.*)"', attr_str.strip())
        if m:
            attr_dict[m.group(1)] = m.group(2)
    return attr_dict


def make_attr_field(attr_dict, attr_list):
    field = list()
    for attr in attr_list:
        if attr in attr_dict:
            field.append(f'{attr} "{attr_dict[attr]}"')
        else:
            field.append(f'{attr} ""')
    return "; ".join(field)


def main():
    if len(sys.argv) != 4:
        print("usage:", sys.argv[0], "cons.gtf tmap.tsv ref.gtf", file=sys.stderr)
        sys.exit(1)

    new_trans_id_dict = dict()
    new_gene_id_dict = dict()
    with open(sys.argv[2]) as fin:
        n_new_trans_dict = defaultdict(int)
        reader = csv.reader(fin, delimiter="\t")
        next(reader)
        for row in reader:
            ref_gene_id = row[0]
            ref_trans_id = row[1]
            class_code = row[2]
            qry_gene_id = row[3]
            qry_trans_id = row[4]

            new_trans_id = qry_trans_id
            if class_code == "=":
                new_trans_id = ref_trans_id
            elif class_code == "~":
                n_new_trans_dict[ref_trans_id] += 1
                new_trans_id = f"{ref_trans_id}.nic.{n_new_trans_dict[ref_trans_id]}"
            elif ref_trans_id != "-":
                n_new_trans_dict[ref_trans_id] += 1
                new_trans_id = f"{ref_trans_id}.nnic.{n_new_trans_dict[ref_trans_id]}"

            new_gene_id = qry_gene_id
            if ref_gene_id != "-":
                new_gene_id = ref_gene_id

            new_trans_id_dict[qry_trans_id] = new_trans_id
            new_gene_id_dict[qry_trans_id] = new_gene_id

    gene_attr_dict = dict()
    with open(sys.argv[3]) as fin:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if f[2] == "gene":
                attr = parse_attr(f[8])
                gene_attr_dict[attr["gene_id"]] = attr

    new_gene_info_dict = dict()
    gene_id_list = list()
    out_dict = dict()
    with open(sys.argv[1]) as fin:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if f[2] != "transcript":
                continue
            attr = parse_attr(f[8])
            new_gene_id = new_gene_id_dict[attr["transcript_id"]]
            start = int(f[3])
            end = int(f[4])
            strand = f[6]
            if new_gene_id in new_gene_info_dict:
                new_gene_info = new_gene_info_dict[new_gene_id]
                if new_gene_info.start > start:
                    new_gene_info.start = start
                if new_gene_info.end < end:
                    new_gene_info.end = end
                if new_gene_info.strand == ".":
                    new_gene_info.strand = strand
                continue

            new_gene_info_dict[new_gene_id] = GeneInfo(f[0], start, end, strand)

            new_attr = {
                "transcript_id": new_trans_id_dict[attr["transcript_id"]],
                "gene_id": new_gene_id_dict[attr["transcript_id"]]
            }
            new_gene_attr = dict()
            if new_attr["gene_id"] in gene_attr_dict:
                new_gene_attr = gene_attr_dict[new_attr["gene_id"]]
                new_gene_attr_field = make_attr_field(gene_attr_dict[new_attr["gene_id"]], ["gene_id", "gene_type", "gene_name"])
            else:
                new_gene_attr = {
                    "gene_id": new_attr["gene_id"],
                    "gene_type": "novel",
                    "gene_name": new_attr["gene_id"],
                }
            new_gene_attr_field = make_attr_field(new_gene_attr, ["gene_id", "gene_type", "gene_name"])
            new_gene_info = new_gene_info_dict[new_gene_attr["gene_id"]]
            out_line = ("\t".join([
                new_gene_info.seq_id,
                "gffcompare",
                "gene",
                str(new_gene_info.start),
                str(new_gene_info.end),
                ".",
                new_gene_info.strand,
                "."
            ]) + "\t" + new_gene_attr_field)
            out_dict[new_gene_attr["gene_id"]] = [out_line]
            gene_id_list.append(new_gene_attr["gene_id"])

    with open(sys.argv[1]) as fin:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if f[2] == "transcript":
                attr = parse_attr(f[8])
                new_attr = {
                    "transcript_id": new_trans_id_dict[attr["transcript_id"]],
                    "gene_id": new_gene_id_dict[attr["transcript_id"]]
                }
                new_attr_field = make_attr_field(new_attr, ["transcript_id", "gene_id"])
            elif f[2] == "exon":
                attr = parse_attr(f[8])
                new_attr = {
                    "transcript_id": new_trans_id_dict[attr["transcript_id"]],
                    "gene_id": new_gene_id_dict[attr["transcript_id"]]
                }
                new_attr_field = make_attr_field(new_attr, ["transcript_id", "gene_id"])
            else:
                continue

            f[1] = "gffcompare"
            out_line = "\t".join(f[0:8]) + "\t" + new_attr_field
            out_dict[new_attr["gene_id"]].append(out_line)

    for gene_id in gene_id_list:
        print("\n".join(out_dict[gene_id]))


if __name__ == "__main__":
    main()
