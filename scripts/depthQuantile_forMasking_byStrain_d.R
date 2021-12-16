args <- commandArgs(TRUE)
strainName <- args[1]
strainFile <- read.table(paste(strainName, ".bedgraph", sep=""), header=FALSE, col.names = c("chrom", "chromstart", "value"))
strainOut <- data.frame(strain=NA, q10=NA, q99=NA)
strainOut$strain <- strainName
strainOut$q10<- quantile(strainFile$value, 0.1)
strainOut$q99<- quantile(strainFile$value, 0.99)
write.table(strainOut, file=paste(strainName, ".cov", sep=""), row.names=F, sep = "\t")
