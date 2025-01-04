#!/usr/bin/python

import logging
import os
import re
import shutil
import sys
from collections import defaultdict
from isorefiner.common import get_version, run_command, func_with_log

logger = logging.getLogger(__name__)

        
def main(args):
    try:
        # Set variables
        raw_reads_files = [os.path.abspath(_) for _ in args.reads]
        raw_genome_file = os.path.abspath(args.genome)
        out_gtf = os.path.abspath(args.out_gtf)
        n_thread = args.threads
        max_mem = args.max_mem
        rnabloom_option = args.rnabloom_option
        gmap_min_cov = args.gmap_min_cov
        gmap_min_idt = args.gmap_min_idt
        gmap_max_intron = args.gmap_max_intron
        gmap_option = args.gmap_option

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
        genome_file = "genome.fasta"
        if os.path.lexists(genome_file):
            genome_file = f"genome_{os.getpid()}.fasta"
        os.symlink(raw_genome_file, genome_file)

        reads_files = list()
        for i, raw_reads_file in enumerate(raw_reads_files, start=1):
            reads_ext = os.path.basename(raw_reads_file).split(".", 1)[1]
            reads_file = f"reads_{i}.{reads_ext}"
            if os.path.lexists(reads_file):
                reads_file = f"reads_{i}_{os.getpid()}.{reads_ext}"
            reads_files.append(reads_file)
            os.symlink(raw_reads_file, reads_file)

        # Main process
        logger.info(f"isorefiner version: {get_version()}")
        logger.info(f"Starting isorefiner {args.command}")
        ## RNA-Bloom
        os.environ["JAVA_TOOL_OPTIONS"] = f"-Xmx{max_mem}"
        run_command(f"rnabloom -t {n_thread} -outdir rnabloom_out {rnabloom_option} -long {' '.join(reads_files)}")
        ## GMAP
        run_command(
            f"gmap_build -D . -d genome_gmap {genome_file}",
            stdout="gmap_build.stdout",
            stderr="gmap_build.stderr"
        )
        run_command(
            " ".join([
                "gmap",
                "-D .",
                "-d genome_gmap",
                "-f 2",
                f"-t {n_thread}",
                f"--min-trimmed-coverage={gmap_min_cov}",
                f"--min-identity={gmap_min_idt}",
                f"--max-intronlength-middle={gmap_max_intron}",
                gmap_option,
                "rnabloom_out/rnabloom.transcripts.fa"
            ]),
            stdout="gmap.gff3",
            stderr="gmap.stderr"
        )
        run_command("gffread -E -T -o gmap.gtf gmap.gff3")
        gtf_filter_long_intron("gmap.gtf", out_gtf, gmap_max_intron)
        shutil.rmtree("genome_gmap")
        logger.info(f"Finished isorefiner {args.command}")

    except Exception as e:
        msg = f"Error: {e}, Exception class: {type(e).__name__}"
        print(msg, file=sys.stderr)
        logger.error(msg)
        sys.exit(1)


@func_with_log
def gtf_filter_long_intron(in_gtf, out_gtf, max_intron_len):
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
            exon_pos_dict[transcript_id].add(f[0], int(f[3]) - 1, int(f[4]))

    rm_set = set()
    for transcript_id, exon_pos in exon_pos_dict.items():
        splice_site_list = exon_pos.get_splice_sites()
        if len(splice_site_list) < 2:
            continue
        for i in range(0, len(splice_site_list), 2):
            if splice_site_list[i + 1] - splice_site_list[i] > max_intron_len:
                rm_set.add(transcript_id)

    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for ln in fin:
            if len(ln) == 0 or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            m = transcript_re.search(f[8])
            if not m:
                continue
            transcript_id = m.group(1)
            if transcript_id not in rm_set:
                print(ln, end="", file=fout)