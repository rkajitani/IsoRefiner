#!/usr/bin/python

import logging
import os
import re
import subprocess
import sys
from collections import defaultdict
from operator import itemgetter
import polars as pl
from isorefiner.common import func_with_log, run_command, filter_bam

logger = logging.getLogger(__name__)


def main(args):
    try:
        # Set variables
        raw_input_gtfs = [os.path.abspath(_) for _ in args.input_gtf]
        raw_genome_file = os.path.abspath(args.genome)
        raw_ref_gtf = os.path.abspath(args.ref_gtf)
        raw_reads_files = [os.path.abspath(_) for _ in args.reads]
        out_gtf = os.path.abspath(args.out_gtf)
        n_thread = args.threads
        max_indel = args.max_indel
        max_clip = args.max_clip
        min_idt = args.min_idt

        # Create and move to working directory
        work_dir = args.work_dir
        os.makedirs(work_dir, exist_ok=True)
        os.chdir(work_dir)

        # Set logger
        logging.basicConfig(
            filename="log.txt",
            level=logging.INFO,
            format="%(asctime)s - PID=%(process)d - %(levelname)s - %(message)s"
        )

        # Set input files and create symbolic links
        input_gtfs = list()
        for i, raw_input_gtf in enumerate(raw_input_gtfs, start=1):
            input_gtf = f"raw_{i}.gtf"
            if os.path.lexists(input_gtf):
                input_gtf = f"raw_{i}_{os.getpid()}.gtf"
            input_gtfs.append(input_gtf)
            os.symlink(raw_input_gtf, input_gtf)

        genome_file = "genome.fasta"
        if os.path.lexists(genome_file):
            genome_file = f"genome_{os.getpid()}.fasta"
        os.symlink(raw_genome_file, genome_file)

        ref_gtf = "ref.gtf"
        if os.path.lexists(ref_gtf):
            ref_gtf = f"ref_{os.getpid()}.gtf"
        os.symlink(raw_ref_gtf, ref_gtf)

        reads_files = list()
        for i, raw_reads_file in enumerate(raw_reads_files, start=1):
            reads_ext = os.path.basename(raw_reads_file).split(".", 1)[1]
            reads_file = f"reads_{i}.{reads_ext}"
            if os.path.lexists(reads_file):
                reads_file = f"reads_{i}_{os.getpid()}.{reads_ext}"
            reads_files.append(reads_file)
            os.symlink(raw_reads_file, reads_file)

        # Main process start
        logger.info(f"Starting isorefiner refine")

        ## Merge step
        run_command(
            f"gffcompare -o merge -p cons -r {ref_gtf} {ref_gtf} {' '.join(input_gtfs)}",
            stdout="gffcompare_merge.stdout",
            stderr="gffcompare_merge.stderr"
        )

        ## Strand-correction step
        run_command(
            f"gffcompare -r {ref_gtf} merge.combined.gtf -o flip",
            stdout="gffcompare_flip.stdout",
            stderr="gffcompare_flip.stderr"
        )
        gtf_flip_strand("merge.combined.gtf", "flip.merge.combined.gtf.tmap", "flipped.gtf")

        ## Correction step based on intron distance 
        run_command(f"gffread -w asm.fa -g {genome_file} flipped.gtf")
        run_command(
            f"minimap2 -ax map-ont --secondary=no -t {n_thread} asm.fa {' '.join(reads_files)} | samtools view -b -F 2308 -",
            stdout="raw.bam"
        )
        run_command(
            f"samtools sort -@ {n_thread} -m 2G raw.bam",
            stdout="sorted.bam",
            stderr="samtools_sort.stderr"
        )
        filter_bam("sorted.bam", "filt.bam", max_indel, max_clip, min_idt)
        run_command(f"samtools index filt.bam")
        run_command(f"samtools coverage filt.bam", stdout="cov_raw.tsv")
        reform_coverage_table("cov_raw.tsv", "cov.tsv")
        gtf_to_exon_bed("flipped.gtf", "exon.bed")
        bed_part_stats("exon.bed", "exon_n_len.tsv")
        gtf_to_intron_bed("flipped.gtf", "intron.bed")
        bed_part_stats("intron.bed", "intron_n_len.tsv")
        run_command(
            f"bedtools intersect -s -wao -a intron.bed -b intron.bed | perl -ane 'print if ($F[3] ne $F[9])'",
            stdout="intron_self_ovl.tsv",
            stderr="bedtools_intersect.stderr"
        )
        intron_pos_dist("intron_self_ovl.tsv", "intron_n_len.tsv", "intron_pos_dist.tsv")
        run_command(
            "|".join([
                f"bedtools getfasta -nameOnly -s -bed intron.bed -fi {genome_file}",
                "perl -pne 's/\([+-]\)//'",
                "seqkit fx2tab",
                'perl -ane \'print(join("\t", ($F[0], substr($F[1], 0, 2), substr($F[1], -2, 2))), "\n")\'',
            ]),
            stdout="intron_edge.tsv",
            stderr="bedtools_getfasta.stderr"
        )
        judge_canonical_intron("intron_edge.tsv", "non_canonical_flag.tsv")
        join_intron_and_dov_tsv("intron_pos_dist.tsv", "cov.tsv", "non_canonical_flag.tsv", "exon_n_len.tsv", "intron_joined.tsv")


        logger.info(f"Finished isorefiner refine")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)


@func_with_log
def gtf_flip_strand(in_gtf, tmap_file, out_gtf):
    target_set = set()
    with open(tmap_file) as fin:
        fin.readline()
        for ln in fin:
            f = ln.rstrip("\n").split("\t")
            if f[2] == "s":
                target_set.add(f[4])

    transcript_re = re.compile(r'transcript_id "([^"]*)"')
    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            m = transcript_re.search(f[8])
            if m and (m.group(1) in target_set):
                if f[6] == "+":
                    f[6] = "-"
                elif f[6] == "-":
                    f[6] = "+"
                print("\t".join(f), file=fout)
            else:
                print(ln, end="", file=fout)


@func_with_log
def reform_coverage_table(in_file, out_file):
    df = (pl.read_csv(in_file, separator="\t")
        .select("#rname", "coverage", "meandepth")
        .rename({"#rname": "qry_id"})
    )
    df.write_csv(out_file, separator="\t")


@func_with_log
def gtf_to_exon_bed(in_gtf, out_bed):
    transcript_re = re.compile(r'transcript_id "([^"]*)"')
    n_exon_dict = defaultdict(int)
    with open(in_gtf) as fin, open(out_bed, "w") as fout:
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
            print(
                f[0],
                int(f[3]) - 1,
                f[4],
                f"{transcript_id},{n_exon_dict[transcript_id]}",
                ".",
                f[6],
                sep="\t",
                file=fout
            )


@func_with_log
def bed_part_stats(in_bed, out_tsv):
    class Transcript:
        def __init__(self):
            self.n = 0
            self.length = 0
        def add(self, length):
            self.n += 1
            self.length += length

    part_id_re = re.compile(r'^(\S+),\d+$')
    part_stat_dict = defaultdict(Transcript)
    with open(in_bed) as fin:
        for ln in fin:
            f = ln.rstrip("\n").split("\t")
            m = part_id_re.match(f[3])
            if not m:
                continue
            transcript_id = m.group(1)
            part_stat_dict[transcript_id].add(int(f[2]) - int(f[1]))

    with open(out_tsv, "w") as fout:
        print("transcript_id", "n_parts", "sum_len", sep="\t", file=fout)
        for transcript_id, part_stat in part_stat_dict.items():
            print(transcript_id, part_stat.n, part_stat.length, sep="\t", file=fout)


@func_with_log
def gtf_to_intron_bed(in_gtf, out_bed):
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
    with open(in_gtf) as fin:
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

    with open(out_bed, "w") as fout:
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
                    sep="\t",
                    file=fout
                )


@func_with_log
def intron_pos_dist(intron_self_ovl_tsv, intron_n_len_tsv, out_tsv):
    class Overlap:
        def __init__(self):
            self.n = 0
            self.dist = 0
        def calc_dist(self, start1, end1, start2, end2):
            self.dist += abs(start1 - start2) + abs(end1 - end2)
            self.n += 1

    n_intron_dict = dict()
    with open(intron_n_len_tsv) as fin:
        fin.readline()
        for ln in fin:
            f = ln.rstrip("\n").split("\t")
            n_intron_dict[f[0]] = int(f[1])

    exon_id_re = re.compile(r'^(\S+),\d+$')
    ovl_dict = defaultdict(lambda: defaultdict(Overlap))
    with open(intron_self_ovl_tsv) as fin:
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

    with open(out_tsv, "w") as fout:
        print(
            "transcript_id",
            "counter_id",
            "intron_pos_dist",
            "n_intron",
            "counter_n_intron",
            sep="\t",
            file=fout
        )
        for transcript_id, each_ovl_dict in ovl_dict.items():
            n_intron = n_intron_dict[transcript_id]
            min_dist_id = ""
            min_dist = 0
            for counter_id, ovl in each_ovl_dict.items():
                if n_intron == ovl.n and (min_dist_id == "" or min_dist > ovl.dist):
                    min_dist = ovl.dist
                    min_dist_id = counter_id
            if min_dist_id != "":
                print(
                    transcript_id,
                    min_dist_id,
                    min_dist,
                    n_intron,
                    n_intron_dict[counter_id],
                    sep="\t",
                    file=fout
                )


@func_with_log
def judge_canonical_intron(in_tsv, out_tsv):
    canonical_set = set([("GT", "AG")])

    exon_id_re = re.compile(r'^(\S+),\d+$')
    id_set = set()
    non_canonical_id_set = set()
    with open(in_tsv) as fin:
        for ln in fin:
            f = ln.rstrip("\n").split("\t")
            m = exon_id_re.match(f[0])
            if not m:
                continue
            transcript_id = m.group(1)

            id_set.add(transcript_id)
            if (f[1], f[2]) not in canonical_set:
                non_canonical_id_set.add(transcript_id)

    with open(out_tsv, "w") as fout:
        print("transcript_id", "non_canonical_flag", sep="\t", file=fout)
        for transcript_id in sorted(id_set):
            non_canonical_flag = 1
            if transcript_id not in non_canonical_id_set:
                non_canonical_flag = 0
            print(transcript_id, non_canonical_flag, sep="\t", file=fout)


@func_with_log
def join_intron_and_dov_tsv(intron_pos_dist_tsv, cov_tsv, non_canonical_flag_tsv, exon_n_len_tsv, out_tsv):
    intron_pos_dist_df = pl.read_csv(intron_pos_dist_tsv, separator="\t")
    cov_df = pl.read_csv(cov_tsv, separator="\t").select("qry_id", "meandepth")
    non_canonical_flag_df = pl.read_csv(non_canonical_flag_tsv, separator="\t")
    exon_n_len_df = pl.read_csv(exon_n_len_tsv, separator="\t").select("transcript_id", "sum_len")

    joined_df = (intron_pos_dist_df
        .join(cov_df, left_on="transcript_id", right_on="qry_id", how="left")
        .rename({"meandepth": "transcript_depth"})
        .join(cov_df, left_on="counter_id", right_on="qry_id", how="left")
        .rename({"meandepth": "counter_depth"})
        .join(non_canonical_flag_df, left_on="counter_id", right_on="transcript_id", how="left")
        .rename({"non_canonical_flag": "counter_non_canonical_flag"})
        .join(non_canonical_flag_df, on="transcript_id", how="left")
        .join(exon_n_len_df, left_on="counter_id", right_on="transcript_id", how="left")
        .rename({"sum_len": "counter_len"})
        .join(exon_n_len_df, on="transcript_id", how="left")
        .rename({"sum_len": "transcript_len"})
        .select(
            "transcript_id",
            "counter_id",
            "intron_pos_dist",
            "n_intron",
            "counter_n_intron",
            "transcript_depth",
            "counter_depth",
            "non_canonical_flag",
            "counter_non_canonical_flag",
            "transcript_len",
            "counter_len"
        )
    )
    joined_df.write_csv(out_tsv, separator="\t")
