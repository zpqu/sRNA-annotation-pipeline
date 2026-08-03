##Date: 2026
##Author: Zhipeng
## Single source of truth for the reference genome AND the annotation strategy
## used by the pipeline. Every R script sources this file (via a
## project-root-relative path) and derives object names, DB paths, the
## sense-annotation rule and output locations from it. The pipeline is
## genome-agnostic: to use a different reference, point `genome` at the
## assembly ID of an existing feature DB and re-run the pipeline.
##
##   genome  : assembly ID used to locate the feature DB and raw files
##             (e.g. "mm39", "hg38"). The feature DB must have been built
##             into DB/rdata_<genome>/ by scripts/R/00_build_DB/0_build_annotation_DB.R.
##   strategy: sense-annotation rule for the small-RNA annotation step.
##             One of
##               "fully-contained"  (default) read fully contained in the feature
##               "union"                     read within feature OR feature within read
##               "any"                       any overlap (>= 1 bp) between read and feature
##               "comparison"                runs fully-contained, union and any and then
##                                           performs the overlap-rule comparison analysis
##             Overridden by the environment variable SMALLRNA_STRATEGY (or, inside a
##             comparison run, by SMALLRNA_SUBSTRATEGY, which the runner sets to the
##             strategy currently being executed).
##
##   samples : sample labels, in order, used by the figure scripts.
##             The annotation steps (1-2) detect samples from the BAM/RData
##             file names automatically; only the figure/comparison scripts
##             use this list.
##
## Location: config/genome.R at the project root. Scripts run from
## scripts/R/<stage>/ and source it as source("../../../config/genome.R").

genome  <- Sys.getenv("SMALLRNA_GENOME", unset = "mm39")  # env var overrides the config default

## project root, resolved regardless of the working directory the script is run from
project.root <- if (basename(getwd()) %in% c("00_build_DB", "01_preprocess",
                                              "02_annotation", "03_figures", "04_summary")) {
  "../../.."
} else if (file.exists("config/genome.R")) {
  "."
} else {
  stop("cannot locate project root from working directory ", getwd())
}

## ---- annotation strategy -----------------------------------------------------
canon.strategy = function(x){
  x = tolower(trimws(x))
  x = gsub("[ _]+", "-", x)
  if (grepl("union", x)) return("union")
  if (grepl("any", x))   return("any")
  if (grepl("compar", x)) return("comparison")
  "fully-contained"
}
strategy.mode <- canon.strategy(Sys.getenv("SMALLRNA_STRATEGY", unset = "fully-contained"))
sub.strategy  <- Sys.getenv("SMALLRNA_SUBSTRATEGY", unset = NA_character_)
is.comparison <- strategy.mode == "comparison"
cur.strategy  <- if (!is.na(sub.strategy)) canon.strategy(sub.strategy) else strategy.mode

## directories (relative to project.root)
db.dir    <- file.path(project.root, "DB", paste0("rdata_", genome))        # built feature GRanges
orig.dir  <- file.path(project.root, "DB", paste0("original_data_", genome)) # raw feature files for the DB build

## output layout:
##   single-strategy run (fully-contained | union | any) -> output/
##   comparison run                                       -> output/comparison/ with one
##      subfolder per strategy (fully_contained/, union/, any/) plus shared step-1
##      results and the comparison-analysis tables/figures at the comparison root.
out.base    <- file.path(project.root, if (is.comparison) "output/comparison" else "output")
strategy.dir <- c("fully-contained" = "fully_contained", "union" = "union", "any" = "any")
out.dir     <- if (is.comparison && !is.na(sub.strategy)) {
                 file.path(out.base, unname(strategy.dir[cur.strategy]))
               } else {
                 out.base
               }

## sample labels (order matters for the figure scripts)
samples <- c("Cumulus-cells", "Granulosa-cells")

## chromosome style enforced on feature seqnames by the DB build ("chr" or "")
chr.style <- "chr"

stopifnot(dir.exists(file.path(project.root, "DB")))
