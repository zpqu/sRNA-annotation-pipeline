##Date: 2026
##Author: Zhipeng
## Single source of truth for the reference genome used by the pipeline.
## Every R script sources this file (via a project-root-relative path) and
## derives object names, DB paths and output locations from it. The pipeline
## is genome-agnostic: to use a different reference, point `genome` at the
## assembly ID of an existing feature DB and re-run the pipeline.
##
##   genome : assembly ID used to locate the feature DB and raw files
##            (e.g. "mm39", "hg38"). The feature DB must have been built
##            into DB/rdata_<genome>/ by scripts/R/00_build_DB/0_build_annotation_DB.R.
##   samples: sample labels, in order, used by the figure scripts.
##            The annotation steps (1-2) detect samples from the BAM/RData
##            file names automatically; only the figure/comparison scripts
##            use this list.
##
## Location: config/genome.R at the project root. Scripts run from
## scripts/R/<stage>/ and source it as source("../../../config/genome.R").

genome  <- Sys.getenv("SMALLRNA_GENOME", unset = "mm39")  # env var overrides the config default

## project root, resolved regardless of the working directory the script is run from
project.root <- if (basename(getwd()) %in% c("00_build_DB", "01_preprocess",
                                              "02_annotation", "03_figures")) {
  "../../.."
} else if (file.exists("config/genome.R")) {
  "."
} else {
  stop("cannot locate project root from working directory ", getwd())
}

## directories (relative to project.root)
db.dir    <- file.path(project.root, "DB", paste0("rdata_", genome))        # built feature GRanges
orig.dir  <- file.path(project.root, "DB", paste0("original_data_", genome)) # raw feature files for the DB build
out.dir   <- file.path(project.root, "output_matmiRNA")                     # main-strategy annotation output
out.dir.pri <- file.path(project.root, "output")                            # primary-miRNA strategy output

## sample labels (order matters for the figure scripts)
samples <- c("Cumulus-cells", "Granulosa-cells")

## chromosome style enforced on feature seqnames by the DB build ("chr" or "")
chr.style <- "chr"

stopifnot(dir.exists(file.path(project.root, "DB")))
