#!/usr/bin/Rscript --vanilla
# Author: Peris, IATA-CSIC

library(ggplot2)
library(gridExtra)
library(dplyr)
library(optparse)

option_list = list(
make_option(c("-i", "--dirFolder"), action="store", type="character", default=NULL,
        help="Direct folder where qaf.tab files are stored",metavar="character"),
make_option(c("-o", "--outputFile"), action="store", type="character", default=NULL,
        help="Name of the output file",metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

data_HTZ_distr <- list.files(opt$dirFolder, pattern = "*.qaf.tab")
lista_p = list()

pdf(paste0(opt$outputFile,"_HTZ_multiplot.pdf"),height=6,width=14)
for(i in 1:length(data_HTZ_distr)){
  f_temp <- read.table(paste0(opt$dirFolder,data_HTZ_distr[i]),sep= "\t", header = T)
  f_temp <- f_temp[nchar(as.character(f_temp$REF))==1,]
  f_temp <- f_temp[nchar(as.character(f_temp$ALT))==1,]
  f_temp <- mutate(f_temp, AFATOT = AFA*TOTC, AFRTOT=AFR*TOTC)
  lista_p[[i]] <- ggplot(f_temp, aes_string(x="AFATOT",y="AFRTOT")) +
      geom_point(size = 1, color = "red") +
      theme_classic() +
      theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), axis.title.x = element_blank(),axis.title.y = element_blank()) +
      xlab("Alternative") + scale_x_continuous(limits = c(0,400)) +
      ylab("Reference") + scale_y_continuous(limits = c(0,400)) + 
      ggtitle(unlist(strsplit(as.character(data_HTZ_distr[i]),".qaf"))[1])
                   
  print(ggplot(f_temp, aes_string(x="AFATOT",y="AFRTOT")) +
      geom_point(size = 1, color = "red") +
      theme_classic() +
      theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), axis.title.x = element_blank(),axis.title.y = element_blank()) +
      xlab("Alternative") + scale_x_continuous(limits = c(0,400)) +
      ylab("Reference") + scale_y_continuous(limits = c(0,400)) + 
      ggtitle(unlist(strsplit(as.character(data_HTZ_distr[i]),".qaf"))[1]))
  names_plots <- c(lista_p, paste0("p",i))
}
print(do.call(grid.arrange,c(lista_p, nrow = 2)))
dev.off()
