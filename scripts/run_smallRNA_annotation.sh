echo -ne "The whole pipeline starts at: "
date

cd bashScripts
./folder_clean.sh


cd ../rScripts
##step 1
echo -ne "Now is formatting the bam files ... at "
date
Rscript 1_format_bam.R

##step 2
echo -ne "Now is annotating small RNAs ... at "
date
Rscript 2_annotation_smallRNA_RMlast.R

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
