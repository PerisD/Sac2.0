#!/usr/bin/Rscript --vanilla
# Author: Miguel Morard and modified by Peris

library(ggplot2)
library(gridExtra)

args <- commandArgs(TRUE)
inputName <- args[1] #PATH of the converted file
chromosomeInfo <- args[2] #"pairwiseDivFile_chromosome.txt"

f = read.table(inputName,sep= "\t", header = T) #Lee la tabla 000.vcf.qaf.tab

chromosome_properties = read.csv(chromosomeInfo,header=TRUE,sep='\t')
allchr <- c(as.character(chromosome_properties$NAME))
chrlen <- list()
for (i in 1:length(chromosome_properties[,1])){
  chrlen[[i]] <- subset(chromosome_properties, chromosome_properties$NAME == allchr[i])
}

lista_p = list() #lista para almacenar los ggplots
subs = list() #lista para almacenar los datos de cada cromosoma
for (i in 1:length(allchr)){
  subs[[i]] <- subset(f, f$CHR == allchr[i])
}

color_species <- c("#ff0000","#cc6600","#008000","#00ffff","#0070c0","#7030a0","#ff00ff")
pdf(paste0(inputName,"_distr.pdf"),height=6,width=14)

print(ggplot(f, aes(x=AFA)) +
  geom_freqpoly(aes(y=..count..),binwidth = 0.03) +
  geom_histogram(binwidth = 0.03) +
  theme_classic() +
  scale_colour_manual(values = mycolor) +
  scale_linetype_manual(values = myshapes) +
  #ggtitle("Total counts") +
  #scale_x_continuous(limits = c(0.0,1.5)) + 
  theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black")))

lista_p = list()
chromCHECK <- c(1,17,33,49,65,81,97)
names_plots <- c()
count <- 1
for(i in 1:length(chrlen)){
  if(i %in% chromCHECK){
    color2plot <- color_species[count]
    lista_p[[i]] <- ggplot(subs[[i]], aes_string(x="POS",y="AFA")) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),plot.title = element_text(hjust = 0.5, colour ="#429221",size = 11),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y = element_blank(),axis.text.y = element_text(face = "bold", size = 16)) +
    geom_point(colour = color2plot, size=1) +
    scale_y_continuous(limits = c(0,1),breaks = c(0,0.5,1)) +
    scale_x_continuous(limits = c(0,chrlen[[i]]$LEN)) +
    xlab(unique(as.character(subs[[i]]$CHR)))
    count <- count + 1
  }else{
    lista_p[[i]] <- ggplot(subs[[i]], mapping = aes_string(x="POS",y="AFA")) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),plot.title = element_text(hjust = 0.5, colour ="#429221",size = 11),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank()) + 
      geom_point(colour = color2plot, size=1) +
      scale_y_continuous(limits = c(0,1)) +
      scale_x_continuous(limits = c(0,chrlen[[i]]$LEN)) +
      xlab(unique(as.character(subs[[i]]$CHR)))
  }
  names_plots <- c(lista_p, paste0("p",i))
}
print(do.call(grid.arrange,c(lista_p, ncol=16)))
#Show a allele frequency distribution (f)

dev.off()

