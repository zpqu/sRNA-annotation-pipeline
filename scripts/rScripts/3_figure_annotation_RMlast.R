##Date: 21/07/2015
##Author: Zhipeng
## This script is used to make plot for size distribution for all small RNA data

library(ggplot2)
library(scales)

###for all small RNAs
all.read.df = NULL
files = list.files(path = "../../output/", pattern = ".annotation.count.txt$")

for(i in seq(along = files)){
      file.name = paste("../../output/", files[i], sep = "")
      class.name = sample.name = files[i]
      sample.name = gsub("\\.\\w+\\.annotation\\.count\\.txt", "", sample.name)
      class.name = gsub("\\.annotation\\.count\\.txt", "", class.name)
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
class.list = c("read", "primiRNA", "piRNA", "snoRNA", "tRNA")
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
      if(class.list[j] == "read"){
       single.class.df$Var1 = factor(single.class.df$Var1, 
       			   levels = c("primiRNA", "snoRNA", "piRNA", "tRNA", "RM", "refGene.NM.exon",
			   "refGene.NM.intron", "lincRNA.exon", "AS.primiRNA", "AS.snoRNA", "AS.piRNA", 
			   "AS.tRNA", "AS.RM", "AS.refGene.NM.exon",
			   "AS.refGene.NM.intron", "AS.lincRNA.exon", "other"))
      }else{
       single.class.df$Var1 = factor(single.class.df$Var1, 
       			   levels = c("CDS", "5UTR", "3UTR", "intron", "up1k",
			   "down1k", "RM", "AS.CDS", "AS.5UTR", "AS.3UTR", 
			   "AS.intron", "AS.up1k", "AS.down1k", "AS.RM", "intergenic"))	
      }

      p.sample.class.barplot = ggplot(data = single.class.df, aes(x = Var1, y = Freq)) +
		   geom_bar(stat = "identity") +
		   xlab("") + ylab("Count") +
		   scale_y_continuous(labels = comma) +		   
		   theme_bw() +
		   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		   facet_wrap(~sample, ncol = 4)
      ggsave(p.sample.class.barplot, file = paste0("../../plots/Fig1a.", class.list[j], "_annotation_barplot.pdf"), width = fig.width, height = fig.height)

      p.sample.class.percent.plot = ggplot(data = single.class.df, aes(x = Var1, y = per)) +
		   geom_bar(stat = "identity") +
		   xlab("") + ylab("") +
		   scale_y_continuous(labels = percent) +
		   theme_bw() +
		   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		   facet_wrap(~sample, ncol = 4)
      ggsave(p.sample.class.percent.plot, file = paste0("../../plots/Fig1b.", class.list[j], "_annotation_barplot.percentage.pdf"), width = fig.width, height = fig.height)
}