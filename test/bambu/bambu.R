args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("usage: Rscript bambu.R n_threads genome.fa annot.gtf bam1 [bam2 ...]\n")
	quit("no", status = 1, runLast = FALSE)
}

suppressPackageStartupMessages(library(bambu))

se <- bambu(reads = args[4:length(args)], ncore = args[1], genome = args[2], annotations = args[3])
writeBambuOutput(se, path = "./bambu_out/")
