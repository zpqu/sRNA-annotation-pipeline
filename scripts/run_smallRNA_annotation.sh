echo -ne "The whole pipeline starts at: "
date

## clean previous results (output subfolders and logs)
for d in ../output/rdata ../output/tables ../output/figures logs; do
    if [ "$(ls -A $d)" ]; then
        rm -r ${d}/*
    fi
done

cd ../R/01_preprocess
##step 1
echo -ne "Now is formatting the bam files ... at "
date
Rscript 1_format_bam.R

cd ../02_annotation
##step 2
echo -ne "Now is annotating small RNAs ... at "
date
Rscript 2_annotation_smallRNA_RMlast.R

cd ../03_figures
##step 3
echo -ne "Now is generating annotation plots ... at "
date
Rscript 3_figure_annotation_RMlast.R

##step 4
echo -ne "Now is generating size distribution plots ... at "
date
Rscript 4_figure_size_RMlast.R

##step 5
echo -ne "Now is generating position distribution plots for tRNAs and snoRNAs ... at "
date
Rscript 5_tRNA_snoRNA_position.widnows.R

echo -ne "The whole pipeline ends at: "
date
