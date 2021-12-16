#!/usr/bin/Rscript --vanilla

library(optparse)
library(dplyr)
library(ggplot2)
library(reshape2)

option_list = list(
make_option(c("-i", "--inputFile"), action="store", type="character", default=NULL,
	help="Input file with haplotype and StrainName column",metavar="character"),
make_option(c("-s", "--Strain"), action="store", type="character", default=NULL,
	help="Strain information table, with a column named StrainName",metavar="character"),
make_option(c("-c", "--column"), action="store", type="character", default=NULL,
	help="Column name to be used for generating the frequency nexus file",metavar="character"),
make_option(c("-o", "--outputFile"), action="store", type="character", default=NULL,
	help="Name of the output file",metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$inputFile)
print(opt$Strain)
HaplotypeFile <- read.csv(opt$inputFile)
StrainInfo <- read.csv(opt$Strain)

haplotypes_info_df <- merge(HaplotypeFile, StrainInfo, by = "StrainName", all = TRUE)
haplotypes_info_df <- haplotypes_info_df[!is.na(haplotypes_info_df$Haplotype),]
write.csv(file="haplotype_info.csv",haplotypes_info_df,row.names = FALSE)

freq_table = table(haplotypes_info_df$Haplotype, haplotypes_info_df[[opt$column]])
freq_table <- as.data.frame.matrix(freq_table)
haplotypes <- rownames(freq_table)
freq_table <- cbind(haplotypes,freq_table)

sink(paste0(opt$outputFile,"_2pasteInNexus.txt"))
cat("Begin Traits;\n")
cat(paste0("Dimensions NTraits=",length(colnames(freq_table))-1,";\n"))
cat("Format labels=yes missing=- separator=Comma;\n")
trait_names <- colnames(freq_table)
trait_names[1] <- "TraitLabels"
trait_names <- paste(trait_names,collapse=" ")
trait_names <- c(trait_names,";\n")
trait_names <- paste(trait_names, collapse = "")
cat(trait_names)
cat("Matrix\n")
counter = 1
while(counter <= length(freq_table[,1])){
	haplotype_now <-  as.character(freq_table[counter,1])
	values_now <- as.character(freq_table[counter,2:ncol(freq_table)])
	values_now <- paste(values_now, collapse = ",")
	values_now <- paste0(values_now,"\n")
	values_now <- paste(haplotype_now,values_now, collapse = " ")
	cat(values_now)
	counter <- counter + 1
}
cat(";\n")
cat("End;")
sink()

