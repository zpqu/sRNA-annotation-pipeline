##Date: 21/07/2015
##Author: Zhipeng
## This script is used to make plot for size distribution for all small RNA data.
## Reads the per-sample *.size.count.txt tables from the strategy output
## directory resolved by config/genome.R (step 2).

library(ggplot2)
library(scales)
rm(list = ls())

source("../../../config/genome.R")
dir.tab = file.path(out.dir, "tables")
dir.fig = file.path(out.dir, "figures")
dir.create(dir.fig, recursive = TRUE, showWarnings = FALSE)

###for all small RNAs
all.read.df = NULL
files = list.files(path = dir.tab, pattern = ".size.count.txt$")
if(length(files) == 0) stop("no *.size.count.txt found in ", dir.tab, " (run step 2 first)")
for(i in seq(along = files)){
      file.name = file.path(dir.tab, files[i])
      class.name = sample.name = files[i]
      sample.name = gsub("\\.\\w+\\.size\\.count\\.txt", "", sample.name)
      sample.name = gsub("\\.(bwa|bowtie2)$", "", sample.name)
      class.name = gsub("\\.size\\.count\\.txt", "", class.name)
      class.name = gsub(".+\\.", "", class.name)

      single.read.df = NULL
      print(paste("Now is processing ...", sample.name, Sys.time()))
      single.read.df = read.delim(file.name, header = T)
      single.read.df$sample = rep(sample.name, nrow(single.read.df))
      single.read.df$class = rep(class.name, nrow(single.read.df))

      single.read.df$per = single.read.df$Freq/sum(single.read.df$Freq)
      all.read.df = rbind(all.read.df, single.read.df)
}

###
class.list = c("read", "matmiRNA", "piRNA", "snoRNA", "tRNA")
sample.num = length(unique(all.read.df$sample))
if(sample.num < 4){
              fig.width = sample.num*4
              fig.height = 4
}else{
        fig.height = (as.integer((sample.num-1)/4)+1)*4
        fig.width = 16
}
for(j in 1:length(class.list)){
      single.class.df = all.read.df[all.read.df$class == class.list[j], ]
      if(nrow(single.class.df) == 0) next
      p.sample.class.barplot = ggplot(data = single.class.df, aes(x = Var1, y = Freq)) +
		   geom_bar(stat = "identity") +
		   xlab("") + ylab("Count") +
		   scale_y_continuous(labels = comma) +
		   theme_bw() +
#		   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		   facet_wrap(~sample, ncol = 4)
      ggsave(p.sample.class.barplot, file = file.path(dir.fig, paste0("Fig2a.", class.list[j], "_size_barplot.pdf")), width = fig.width, height = fig.height)

      p.sample.class.percent.plot = ggplot(data = single.class.df, aes(x = Var1, y = per)) +
		   geom_bar(stat = "identity") +
		   xlab("") + ylab("") +
		   scale_y_continuous(labels = percent) +
		   theme_bw() +
#		   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		   facet_wrap(~sample, ncol = 4)
      ggsave(p.sample.class.percent.plot, file = file.path(dir.fig, paste0("Fig2b.", class.list[j], "_size_barplot.percentage.pdf")), width = fig.width, height = fig.height)
}
