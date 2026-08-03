#!/usr/bin/env bash
## Small-RNA annotation pipeline runner.
## Usage: run_smallRNA_annotation.sh [GENOME] [STRATEGY]
##
##   GENOME   reference genome assembly, default mm39 (also SMALLRNA_GENOME).
##            Must match a feature DB built into ../DB/rdata_<GENOME>/
##   STRATEGY sense-annotation rule:
##              fully-contained  (default) read fully contained in the feature
##              union                     read within feature OR feature within read
##              any                       any overlap (>= 1 bp) between read and feature
##              comparison                runs fully-contained, union and any (each
##                                        through steps 2/3/4/5/9) and then performs
##                                        the overlap-rule comparison (step 6) and the
##                                        end-of-pipeline summary (step 10)
##
## Output layout (config/genome.R):
##   single strategy -> output/          (tables/figures/rdata)
##   comparison      -> output/comparison/  (shared step-1 at the root, one
##                                           subfolder per strategy + comparison
##                                           analysis and summary at the root)
set -e
GENOME="${1:-${SMALLRNA_GENOME:-mm39}}"
STRATEGY="${2:-${SMALLRNA_STRATEGY:-fully-contained}}"
case "$STRATEGY" in
  fully-contained|fully_contained|within) STRATEGY="fully-contained" ;;
  union)        STRATEGY="union" ;;
  any)          STRATEGY="any" ;;
  comparison|compar) STRATEGY="comparison" ;;
  *) echo "ERROR: unknown STRATEGY '$STRATEGY' (use fully-contained|union|any|comparison)"; exit 1 ;;
esac
export SMALLRNA_GENOME=$GENOME
export SMALLRNA_STRATEGY=$STRATEGY
if [ ! -d "../DB/rdata_${GENOME}" ]; then
    echo "ERROR: feature DB ../DB/rdata_${GENOME} not found."
    echo "Build it first: cd scripts/R/00_build_DB && Rscript 0_build_annotation_DB.R"
    exit 1
fi
echo -ne "Reference genome: $GENOME, strategy: $STRATEGY. The whole pipeline starts at: "
date

## clean previous results (output base for this run + logs)
BASE=../output
if [ "$STRATEGY" == "comparison" ]; then BASE=../output/comparison; fi
rm -rf "$BASE"
mkdir -p "$BASE/tables" "$BASE/rdata" "$BASE/figures"
rm -rf logs
mkdir -p logs

## run one pipeline script from its stage directory, logging to scripts/logs/
run_R() {
    local stage_dir="$1" script="$2" log="$3"
    echo -ne "$log starts at: "
    date
    ( cd "$stage_dir" && Rscript "$script" > "../../logs/$log" 2>&1 )
    echo -ne "$log finished at: "
    date
}

if [ "$STRATEGY" == "comparison" ]; then
    ## step 1: shared pre-processing (runs once)
    run_R R/01_preprocess 1_format_bam.R step1.log

    ## steps 2/3/4/5/9 for each strategy
    for SUB in fully-contained union any; do
        export SMALLRNA_SUBSTRATEGY=$SUB
        echo "===== strategy: $SUB ====="
        run_R R/02_annotation 2_annotation_smallRNA.R "step2_${SUB}.log"
        run_R R/03_figures 3_figure_annotation.R "step3_${SUB}.log"
        run_R R/03_figures 4_figure_size.R "step4_${SUB}.log"
        run_R R/03_figures 5_tRNA_snoRNA_position.windows.R "step5_${SUB}.log"
        run_R R/03_figures 9_abundance_by_category.R "step9_${SUB}.log"
    done
    unset SMALLRNA_SUBSTRATEGY

    ## step 6: overlap-rule comparison analysis
    run_R R/03_figures 6_compare_overlap_rules.R step6.log
else
    ## step 1
    run_R R/01_preprocess 1_format_bam.R step1.log

    ## step 2
    run_R R/02_annotation 2_annotation_smallRNA.R step2.log

    ## steps 3/4/5
    run_R R/03_figures 3_figure_annotation.R step3.log
    run_R R/03_figures 4_figure_size.R step4.log
    run_R R/03_figures 5_tRNA_snoRNA_position.windows.R step5.log

    ## step 9
    run_R R/03_figures 9_abundance_by_category.R step9.log
fi

## step 10: end-of-pipeline summary
run_R R/04_summary 10_pipeline_summary.R step10.log

echo -ne "The whole pipeline ends at: "
date
