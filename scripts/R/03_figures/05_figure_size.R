##Date: 21/07/2015
##Author: Zhipeng
## This script is used to make plot for size distribution for all small RNA data.
## Reads the consolidated all-reads count table Table_02b from the strategy
## output directory resolved by the shared bootstrap (init.R) and plots, for each
## annotation class, the read-size distribution (faceted by sample).
##
## Outputs (step 05):
##   figures/Figure_05a.<class>_size_barplot.pdf              (counts, all reads)
##   figures/Figure_05b.<class>_size_barplot.percentage.pdf   (percentage)

library(ggplot2)
library(scales)
rm(list = ls())

source("../lib/init.R")
dir.tab = file.path(out.dir, "tables")
dir.fig = file.path(out.dir, "figures")
dir.create(dir.fig, recursive = TRUE, showWarnings = FALSE)

###consolidated all-reads count table from step 02 (Table_02b)
t2b.file = file.path(dir.tab, "Table_02b_annotation_count_all_reads.csv")
if(!file.exists(t2b.file)) stop("no Table_02b_annotation_count_all_reads.csv found in ", dir.tab, " (run step 02 first)")
all.read.df = data.table::as.data.table(read.csv(t2b.file))

class.list = c("read", "matmiRNA", "piRNA", "snoRNA", "tRNA")
size.levels = sort(unique(as.integer(all.read.df$item[all.read.df$category %in%
                        paste0(class.list, ".size")])))

sample.num = length(unique(all.read.df$sample))
if(sample.num < 4){
              fig.width = sample.num*4
              fig.height = 4
}else{
        fig.height = (as.integer((sample.num-1)/4)+1)*4
        fig.width = 16
}
for(j in 1:length(class.list)){
      single.class.df = all.read.df[all.read.df$category == paste0(class.list[j], ".size"), ]
      if(nrow(single.class.df) == 0) next
      single.class.df$per = single.class.df$Freq/sum(single.class.df$Freq)
      single.class.df$item = factor(single.class.df$item, levels = as.character(size.levels))
      p.sample.class.barplot = ggplot(data = single.class.df, aes(x = item, y = Freq)) +
		   geom_bar(stat = "identity") +
		   xlab("Read size (nt)") + ylab("Count") +
		   scale_y_continuous(labels = comma) +
		   theme_bw() +
		   facet_wrap(~sample, ncol = 4)
      ggsave(p.sample.class.barplot, file = file.path(dir.fig, paste0("Figure_05a.", class.list[j], "_size_barplot.pdf")), width = fig.width, height = fig.height)

      p.sample.class.percent.plot = ggplot(data = single.class.df, aes(x = item, y = per)) +
		   geom_bar(stat = "identity") +
		   xlab("Read size (nt)") + ylab("Percentage") +
		   scale_y_continuous(labels = percent) +
		   theme_bw() +
		   facet_wrap(~sample, ncol = 4)
      ggsave(p.sample.class.percent.plot, file = file.path(dir.fig, paste0("Figure_05b.", class.list[j], "_size_barplot.percentage.pdf")), width = fig.width, height = fig.height)
}
