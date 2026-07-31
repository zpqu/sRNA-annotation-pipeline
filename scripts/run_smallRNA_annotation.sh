echo -ne "The whole pipeline starts at: "
date

## clean previous results (output subfolders and logs)
for d in ../output/rdata ../output/tables ../output/figures \
         ../output_matmiRNA/rdata ../output_matmiRNA/tables ../output_matmiRNA/figures logs; do
    if [ "$(ls -A $d)" ]; then
        rm -r ${d}/*
    fi
done

cd R/01_preprocess
##step 1
echo -ne "Step 1 (format bam) starts at: "
date
Rscript 1_format_bam.R > ../../logs/step1.log 2>&1
echo -ne "Step 1 (format bam) finished at: "
date

cd ../02_annotation
##step 2
echo -ne "Step 2 (annotation) starts at: "
date
Rscript 2_annotation_smallRNA_RMlast.R > ../../logs/step2.log 2>&1
echo -ne "Step 2 (annotation) finished at: "
date

##step 2b (mature-miRNA strategy)
echo -ne "Step 2b (mature-miRNA annotation) starts at: "
date
Rscript 2b_annotation_smallRNA_matmiRNA_RMlast.R > ../../logs/step2b.log 2>&1
echo -ne "Step 2b (mature-miRNA annotation) finished at: "
date

cd ../03_figures
##step 3
echo -ne "Step 3 (annotation plots) starts at: "
date
Rscript 3_figure_annotation_RMlast.R > ../../logs/step3.log 2>&1
echo -ne "Step 3 (annotation plots) finished at: "
date

##step 4
echo -ne "Step 4 (size distribution plots) starts at: "
date
Rscript 4_figure_size_RMlast.R > ../../logs/step4.log 2>&1
echo -ne "Step 4 (size distribution plots) finished at: "
date

##step 5
echo -ne "Step 5 (tRNA/snoRNA position plots) starts at: "
date
Rscript 5_tRNA_snoRNA_position.widnows.R > ../../logs/step5.log 2>&1
echo -ne "Step 5 (tRNA/snoRNA position plots) finished at: "
date

##step 6 (mature vs primary miRNA strategy comparison)
echo -ne "Step 6 (miRNA strategy comparison) starts at: "
date
Rscript 6_compare_miRNA_strategies.R > ../../logs/step6.log 2>&1
echo -ne "Step 6 (miRNA strategy comparison) finished at: "
date

echo -ne "The whole pipeline ends at: "
date
