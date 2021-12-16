args <- commandArgs(TRUE)
strainName <- args[1]
strain <- read.table(paste(strainName, sep=""), header=T)
stepSize <- as.numeric(args[2])
outputfolder <- args[3]

strainName <- unlist(strsplit(strainName,"/"))
strainName <- strainName[length(strainName)]
strainName <- tools::file_path_sans_ext(strainName)

completeChromList <- unlist(list(unique(strain$chr)))

genome <- data.frame(chr = NA, start = NA, end = NA)
genome$chr <- "wholeGenome"
genome$start <- 0
genome$end <- strain[length(strain$genomePos), 3]
genome$SNPnumber <- length(strain$Genotype)
genome$SNPprop <- length(strain$Genotype)/strain[length(strain$genomePos), 3]
homozygous <- subset(strain, Genotype==1)
genome$homozygousCount <- length(homozygous$chr)
genome$SNPhomozygousProp <- length(homozygous$chr)/length(strain$Genotype)
genome$totalHomozygousProp <- length(homozygous$chr)/strain[length(strain$genomePos), 3] 
heterozygous <- subset(strain, Genotype==2) 
genome$heterozygousCount <- length(heterozygous$chr)
genome$SNPheterozygousProp <- length(heterozygous$chr)/length(strain$Genotype)
genome$totalHeterozygousProp <- length(heterozygous$chr)/strain[length(strain$genomePos), 3]

output <- data.frame()
for (i in 1:length(completeChromList)) {
  chr_name <- toString(completeChromList[[i]])
  #print(chr_name)
  chr <- subset(strain, chr==chr_name)
  chrSubset <- data.frame(chr = NA, start = NA, end = NA)
  chrSubset$chr <- chr_name
  chrSubset$start <- 0
  chrSubset$end <- chr[length(chr$genomePos), 2]
  chrSubset$SNPnumber <- length(chr$Genotype)
  chrSubset$SNPprop <- length(chr$Genotype)/strain[length(chr$genomePos), 2]
  homozygous <- subset(chr, Genotype==1)
  chrSubset$homozygousCount <- length(homozygous$chr)
  chrSubset$SNPhomozygousProp <- length(homozygous$chr)/length(chr$Genotype)
  chrSubset$totalHomozygousProp <- length(homozygous$chr)/chr[length(chr$genomePos), 2] 
  heterozygous <- subset(chr, Genotype==2) 
  chrSubset$heterozygousCount <- length(heterozygous$chr)
  chrSubset$SNPheterozygousProp <- length(heterozygous$chr)/length(chr$Genotype)
  chrSubset$totalHeterozygousProp <- length(heterozygous$chr)/chr[length(chr$genomePos), 2]
  genome <- rbind(genome, chrSubset)
  
  chr.out <- data.frame()
  counter = 0
  start = 0
  
  if (0 < length(chr[,1])) {
    while (start < chr[length(chr$chr),2]) {
      window.out <- data.frame(chr = NA, start = NA, end = NA)
      window <- subset(chr, (chrPos>start & chrPos<=counter+stepSize))
      window.out$chr <- chr_name  
      window.out$start <- counter
      window.out$end <- counter + stepSize
      window.out$SNPnumber <- length(window$Genotype)
      window.out$SNPprop <- length(window$Genotype)/stepSize
      homozygous <- subset(window, Genotype==1)
      window.out$homozygousCount <- length(homozygous$chr)
      window.out$SNPhomozygousProp <- length(homozygous$chr)/length(window$Genotype)
      window.out$totalHomozygousProp <- length(homozygous$chr)/stepSize 
      heterozygous <- subset(window, Genotype==2)      
      if (0 < length(heterozygous$chr)) {        
        window.out$heterozygousCount <- length(heterozygous$chr)
        window.out$SNPheterozygousProp <- length(heterozygous$chr)/length(window$Genotype)
        window.out$totalHeterozygousProp <- length(heterozygous$chr)/stepSize
      } else {
        window.out$heterozygousCount <- 0
        window.out$SNPheterozygousProp <- 0
        window.out$totalHeterozygousProp <- 0
      }
      start <- counter + stepSize 
      chr.out <- rbind(chr.out, window.out)   
      counter = counter + stepSize
    }   
  } else {
    chr.out$chrom <- chr_name
    chr.out$start <- 0
    chr.out$end <- 0 + stepSize
    chr.out$SNPnumber <-  0
    chr.out$SNPprop <- 0
    chr.out$homozygousCount <- 0
    chr.out$SNPhomozygousProp <- 0
    chr.out$totalHomozygousProp <- 0
    chr.out$heterozygousCount <- 0
    chr.out$SNPheterozygousProp <- 0
    chr.out$totalHeterozygousProp <- 0    
  }  
  output <- rbind(output, chr.out)
}
output <- cbind(Genome_Pos = seq(0, (length(output$chr)-1)*stepSize, stepSize), output)
outputFileName <- paste(paste0(outputfolder,strainName), stepSize, "genotypeSummary.txt", sep="_")
write.table(format(output, scientific=FALSE), file=outputFileName, row.names=F, sep = "\t")

summaryDB <- data.frame(chr = NA, start = NA, end = NA)
summaryDB$chr <- "min"
summaryDB$start <- 0
summaryDB$end <- strain[length(strain$genomePos), 3]
summaryDB$SNPnumber <- min(output$SNPnumber)
summaryDB$SNPprop <- min(output$SNPprop)
summaryDB$homozygousCount <- min(output$homozygousCount)
summaryDB$SNPhomozygousProp <- min(output$SNPhomozygousProp)
summaryDB$totalHomozygousProp <- min(output$totalHomozygousProp)
summaryDB$heterozygousCount <- min(output$heterozygousCount)
summaryDB$SNPheterozygousProp <- min(output$SNPheterozygousProp)
summaryDB$totalHeterozygousProp <- min(output$totalHeterozygousProp)
genome <- rbind(summaryDB, genome)
summaryDB <- data.frame(chr = NA, start = NA, end = NA)
summaryDB$chr <- "mean"
summaryDB$start <- 0
summaryDB$end <- strain[length(strain$genomePos), 3]
summaryDB$SNPnumber <- mean(output$SNPnumber)
summaryDB$SNPprop <- mean(output$SNPprop)
summaryDB$homozygousCount <- mean(output$homozygousCount)
summaryDB$SNPhomozygousProp <- mean(output$SNPhomozygousProp)
summaryDB$totalHomozygousProp <- mean(output$totalHomozygousProp)
summaryDB$heterozygousCount <- mean(output$heterozygousCount)
summaryDB$SNPheterozygousProp <- mean(output$SNPheterozygousProp)
summaryDB$totalHeterozygousProp <- mean(output$totalHeterozygousProp)
genome <- rbind(summaryDB, genome)
summaryDB <- data.frame(chr = NA, start = NA, end = NA)
summaryDB$chr <- "max"
summaryDB$start <- 0
summaryDB$end <- strain[length(strain$genomePos), 3]
summaryDB$SNPnumber <- max(output$SNPnumber)
summaryDB$SNPprop <- max(output$SNPprop)
summaryDB$homozygousCount <- max(output$homozygousCount)
summaryDB$SNPhomozygousProp <- max(output$SNPhomozygousProp)
summaryDB$totalHomozygousProp <- max(output$totalHomozygousProp)
summaryDB$heterozygousCount <- max(output$heterozygousCount)
summaryDB$SNPheterozygousProp <- max(output$SNPheterozygousProp)
summaryDB$totalHeterozygousProp <- max(output$totalHeterozygousProp)
genome <- rbind(summaryDB, genome)
avgFileName <- paste(paste0(outputfolder,strainName), stepSize, "avgGenotype.txt", sep="_")
write.table(genome, file=avgFileName, row.names=F, sep = "\t")

library("ggplot2",lib.loc="/opt/bifxapps/R/library/")
library("grid",lib.loc="/opt/bifxapps/R/library/")
#library(dplyr)
#library("dplyr",lib.loc="/opt/bifxapps/R/library/")
outputPrefix <- paste0(outputfolder,strainName)
window = stepSize

strain <- output
strain[is.na(strain)] <- 0
completeChromList <- unlist(list(unique(strain$chr)))

plot <- ggplot()
xaxis <- scale_x_continuous()
vertLines <- geom_vline()
chrBreaks <- c()
chrLabel <- c()
lineBreaks <- c()

for (i in 1:length(completeChromList)) {
  chr_name = toString(completeChromList[[i]])
  chr_strain <- subset(strain, chr==chr_name)
  chrLabel <- append(chrLabel, chr_name)
  lineBreaks <- append(lineBreaks, chr_strain[1,1])
  subLen <- length(chr_strain[,1]) 
  midpoint <- chr_strain[1,1]+(chr_strain[subLen,3]/2)
  chrBreaks <- append(chrBreaks, midpoint)
}

endPos <- strain[length(strain[,1]),1]
xaxis <- scale_x_continuous(breaks=chrBreaks, labels=chrLabel, name="Genome Position", limits=c(0,endPos))
vertLines <- geom_vline(xintercept = lineBreaks, linetype=2)
line <- geom_abline(intercept=0, slope=0)
midLine <- geom_abline(intercept=50, slope=0, linetype="dotted")

plotSNPnumber <-geom_point(data=strain, aes(x=Genome_Pos, y=SNPnumber,colour="SNP total"))
plotSNPline <- geom_line(data=strain, aes(x=Genome_Pos, y=SNPnumber, colour="SNP total"), size=.75)
plotSNPprop <-geom_point(data=strain, aes(x=Genome_Pos, y=SNPprop*100,colour="SNP total"))
plotHomozygousCount <-geom_point(data=strain, aes(x=Genome_Pos, y=homozygousCount, colour="homozygous"))
plotSNPhomozygousProp <-geom_point(data=strain, aes(x=Genome_Pos, y=SNPhomozygousProp*100, colour="homozygous"))
plotTotalHomozygousProp <-geom_point(data=strain, aes(x=Genome_Pos, y=totalHomozygousProp*100, colour="homozygous"))
plotHeterozygousCount <-geom_point(data=strain, aes(x=Genome_Pos, y=heterozygousCount, colour="heterozygous"))
plotSNPheterozygousProp <-geom_point(data=strain, aes(x=Genome_Pos, y=SNPheterozygousProp*100, colour="heterozygous"))
plotNormSNPheterozygousProp <- geom_ribbon(data=strain, aes(x=Genome_Pos, ymin=0, ymax=SNPheterozygousProp*100, fill="heterozygous"))
plotRibbonSNPheterozygousProp <- geom_ribbon(data=strain, aes(x=Genome_Pos, ymin=0, ymax=SNPheterozygousProp*max(strain$SNPnumber), fill="heterozygous"))
plotTotalHeterozygousProp <-geom_point(data=strain, aes(x=Genome_Pos, y=totalHeterozygousProp*100, colour="heterozygous"))
plotFillSNPtotal <- geom_ribbon(data=strain, aes(x=Genome_Pos, ymin=0, ymax=SNPnumber, fill="SNP total"))
plotFillHetTotal <- geom_ribbon(data=strain, aes(x=Genome_Pos, ymin=0, ymax=heterozygousCount, fill="heterozygous"))
colors <- c("SNP total"= "blue", "homozygous"="red", "heterozygous"="purple")
halfLine <- geom_hline(yintercept=max(strain$SNPnumber)/2, colour="gray", linetype="dotted")

pdf(paste(outputPrefix, "_heterozygosity_", toString(window), ".pdf", sep=""), width=14)

legendLine <- scale_colour_manual(name="Genotype", values=colors)
legendLineFill <- scale_fill_manual(name="Genotype", values=colors)

yaxis <- scale_y_continuous(name="Total Count")
plotTitle <- ggtitle(paste(outputPrefix, "SNP count", sep=" "))
#plot(plot+plotTitle+xaxis+yaxis+line+plotSNPnumber+plotHomozygousCount+plotHeterozygousCount+vertLines+ theme_classic())

yaxis <- scale_y_continuous(name="Precent of window")
plotTitle <- ggtitle(paste(outputPrefix, "precent of window", sep=" "))
#plot(plot+plotTitle+xaxis+yaxis+line+plotSNPprop+plotTotalHomozygousProp+plotTotalHeterozygousProp+vertLines+ theme_classic())

yaxis <- scale_y_continuous(name="Precent of SNPs")
plotTitle <- ggtitle(paste(outputPrefix, "precent of SNPs", sep=" "))
#plot(plot+plotTitle+xaxis+yaxis+line+plotSNPhomozygousProp+plotSNPheterozygousProp+vertLines+midLine+ theme_classic())
yaxis <- scale_y_continuous(name="Precent of SNPs", limits=c(0,100))
#plot(plot+plotTitle+xaxis+yaxis+line+plotNormSNPheterozygousProp+vertLines+midLine+ theme_classic())

yaxis <- scale_y_continuous(name="Total Count")
plotTitle <- ggtitle(paste(outputPrefix, "Proportion of SNPs", sep=" "))
plot(plot+plotTitle+xaxis+yaxis+plotRibbonSNPheterozygousProp+plotSNPline+vertLines+line+halfLine+legendLine+legendLineFill+ theme_classic())

yaxis <- scale_y_continuous(name="Total Count")
plotTitle <- ggtitle(paste(outputPrefix, "SNP totals", sep=" "))
plot(plot+plotTitle+xaxis+yaxis+plotFillSNPtotal+plotFillHetTotal+vertLines+line+halfLine+legendLineFill+ theme_classic())


for (i in 1:length(completeChromList)) {
  chr_name = toString(completeChromList[[i]])
  chrSub <- subset(strain, chr==chr_name)
  plotSNPnumber <-geom_point(data=chrSub, aes(x=start, y=SNPnumber, colour="SNP total"))
  plotSNPline <- geom_line(data=chrSub, aes(x=start, y=SNPnumber, colour="SNP total"), size=.75)
  plotSNPprop <-geom_point(data=chrSub, aes(x=start, y=SNPprop*100, colour="SNP total"))
  plotHomozygousCount <-geom_point(data=chrSub, aes(x=start, y=homozygousCount, colour="homozygous"))
  plotSNPhomozygousProp <-geom_point(data=chrSub, aes(x=start, y=SNPhomozygousProp*100, colour="homozygous"))
  plotTotalHomozygousProp <-geom_point(data=chrSub, aes(x=start, y=totalHomozygousProp*100, colour="homozygous"))
  plotHeterozygousCount <-geom_point(data=chrSub, aes(x=start, y=heterozygousCount, colour="heterozygous"))
  plotSNPheterozygousProp <-geom_point(data=chrSub, aes(x=start, y=SNPheterozygousProp*100, colour="heterozygous"))
  plotNormSNPheterozygousProp <- geom_ribbon(data=chrSub, aes(x=start, ymin=0, ymax=SNPheterozygousProp*100, fill="heterozygous"))
  plotRibbonSNPheterozygousProp <- geom_ribbon(data=chrSub, aes(x=start, ymin=0, ymax=SNPheterozygousProp*max(strain$SNPnumber), fill="heterozygous"))
  plotTotalHeterozygousProp <-geom_point(data=chrSub, aes(x=start, y=totalHeterozygousProp*100, colour="heterozygous")) 
  plotFillSNPtotal <- geom_ribbon(data=chrSub, aes(x=start, ymin=0, ymax=SNPnumber, fill="SNP total"))
  plotFillHetTotal <- geom_ribbon(data=chrSub, aes(x=start, ymin=0, ymax=heterozygousCount, fill="heterozygous"))
  xaxis <- scale_x_continuous(breaks=seq(from=0, to=1200000, by=50000), labels=c("0kb", "50kb", "100kb", "150kb", "200kb", "250kb", "300kb", "350kb", "400kb", "450kb", "500kb", "550kb", "600kb", "650kb", "700kb", "750kb", "800kb", "850kb", "900kb", "950kb", "1Mbp", "1.05Mbp", "1.1Mbp", "1.15Mbp", "1.2Mbp"), name="Chromosome Position")
  yaxis <- scale_y_continuous(name="Precent of window")
  plotTitle <- ggtitle(paste(chr_name, outputPrefix, "precent of window", sep=" "))
  #plot(plot+plotTitle+xaxis+yaxis+line+plotSNPprop+plotTotalHomozygousProp+plotTotalHeterozygousProp+ theme_classic())
  
  yaxis <- scale_y_continuous(name="Precent of SNPs")
  plotTitle <- ggtitle(paste(chr_name, outputPrefix, "precent of SNPs", sep=" "))
  #plot(plot+plotTitle+xaxis+yaxis+line+plotSNPhomozygousProp+plotSNPheterozygousProp+midLine+ theme_classic())
  yaxis <- scale_y_continuous(name="Precent of SNPs", limits=c(0,100))
  #plot(plot+plotTitle+xaxis+yaxis+line+plotNormSNPheterozygousProp+midLine+ theme_classic())
  yaxis <- scale_y_continuous(name="Total Count")
  plotTitle <- ggtitle(paste(chr_name, outputPrefix, "Proportion of SNPs", sep=" "))
  plot(plot+plotTitle+xaxis+yaxis+plotRibbonSNPheterozygousProp+plotSNPline+line+halfLine+legendLine+legendLineFill+ theme_classic())
  yaxis <- scale_y_continuous(name="Total Count")
  plotTitle <- ggtitle(paste(chr_name, outputPrefix, "SNP totals", sep=" "))
  plot(plot+plotTitle+xaxis+yaxis+plotFillSNPtotal+plotFillHetTotal+line+halfLine+legendLineFill+ theme_classic())
}

#Region of interest: chrX 446000 460000
# if (length(completeChromList)<19) {
#   chrX <- subset(strain, chr=="chrX")
# } else {
#   chrX <- subset(strain, chr=="Seub_chrX")
# }
# upstream <- subset(chrX, (start>=432000 & end<=446000))
# region <- subset(chrX, (start>=446000 & end<=460000))
# downstream <- subset(chrX, (start>=460000 & end<=474000))
# 
# plotUpstreamSNPnumber <-geom_point(data=upstream, aes(x=start, y=SNPnumber, colour="flanking SNP total"))
# plotRegionSNPnumber <-geom_point(data=region, aes(x=start, y=SNPnumber, colour="region SNP total"))
# plotDownstreamSNPnumber <-geom_point(data=downstream, aes(x=start, y=SNPnumber, colour="flanking SNP total"))
# plotUpstreamSNPprop <-geom_point(data=upstream, aes(x=start, y=SNPprop*100,colour="flanking SNP total"))
# plotRegionSNPprop <-geom_point(data=region, aes(x=start, y=SNPprop*100, colour="region SNP total"))
# plotDownstreamSNPprop <-geom_point(data=downstream, aes(x=start, y=SNPprop*100, colour="flanking SNP total"))
# plotUpstreamHomozygousCount <-geom_point(data=upstream, aes(x=start, y=homozygousCount, colour="flanking homozygous"))
# plotRegionHomozygousCount <-geom_point(data=region, aes(x=start, y=homozygousCount, colour="region homozygous"))
# plotDownstreamHomozygousCount <-geom_point(data=downstream, aes(x=start, y=homozygousCount, colour="flanking homozygous"))
# plotUpstreamSNPhomozygousProp <-geom_point(data=upstream, aes(x=start, y=SNPhomozygousProp*100, colour="flanking homozygous"))
# plotRegionSNPhomozygousProp <-geom_point(data=region, aes(x=start, y=SNPhomozygousProp*100, colour="region homozygous"))
# plotDownstreamSNPhomozygousProp <-geom_point(data=downstream, aes(x=start, y=SNPhomozygousProp*100, colour="flanking homozygous"))
# plotUpstreamTotalHomozygousProp <-geom_point(data=upstream, aes(x=start, y=totalHomozygousProp*100, colour="flanking homozygous"))
# plotRegionTotalHomozygousProp <-geom_point(data=region, aes(x=start, y=totalHomozygousProp*100, colour="region homozygous"))
# plotDownstreamTotalHomozygousProp <-geom_point(data=downstream, aes(x=start, y=totalHomozygousProp*100, colour="flanking homozygous"))
# plotUpstreamHeterozygousCount <-geom_point(data=upstream, aes(x=start, y=heterozygousCount, colour="flanking heterozygous"))
# plotRegionHeterozygousCount <-geom_point(data=region, aes(x=start, y=heterozygousCount, colour="region heterozygous"))
# plotDownstreamHeterozygousCount <-geom_point(data=downstream, aes(x=start, y=heterozygousCount, colour="flanking heterozygous"))
# plotUpstreamSNPheterozygousProp <-geom_point(data=upstream, aes(x=start, y=SNPheterozygousProp*100, colour="flanking heterozygous"))
# plotRegionSNPheterozygousProp <-geom_point(data=region, aes(x=start, y=SNPheterozygousProp*100, colour="region heterozygous"))
# plotDownstreamSNPheterozygousProp <-geom_point(data=downstream, aes(x=start, y=SNPheterozygousProp*100, colour="flanking heterozygous"))
# plotUpstreamTotalHeterozygousProp <-geom_point(data=upstream, aes(x=start, y=totalHeterozygousProp*100, colour="flanking heterozygous"))
# plotRegionTotalHeterozygousProp <-geom_point(data=region, aes(x=start, y=totalHeterozygousProp*100, colour="region heterozygous"))
# plotDownstreamTotalHeterozygousProp <-geom_point(data=downstream, aes(x=start, y=totalHeterozygousProp*100, colour="flanking heterozygous"))
# 
# vertLines <- geom_vline(xintercept = c(446000,460000))
# xaxis <- scale_x_continuous(name="Seub ChrX", limits=c(432000,474000))
# yaxis <- scale_y_continuous(name="Total Count")
# plotTitle <- ggtitle(paste(outputPrefix, "Saaz vs Frohberg region SNP count", sep=" "))
# plot(plot+plotTitle+xaxis+yaxis+line+plotUpstreamSNPnumber+plotRegionSNPnumber+plotDownstreamSNPnumber+plotUpstreamHomozygousCount+plotRegionHomozygousCount+plotDownstreamHomozygousCount+plotUpstreamHeterozygousCount+plotRegionHeterozygousCount+plotDownstreamHeterozygousCount+vertLines+ theme_classic())
# 
# yaxis <- scale_y_continuous(name="Precent of window")
# plotTitle <- ggtitle(paste(outputPrefix, "Saaz vs Frohberg region precent of window", sep=" "))
# plot(plot+plotTitle+xaxis+yaxis+line+plotUpstreamSNPprop+plotRegionSNPprop+plotDownstreamSNPprop+plotUpstreamTotalHomozygousProp+plotRegionTotalHomozygousProp+plotDownstreamTotalHomozygousProp+plotUpstreamTotalHeterozygousProp+plotRegionTotalHeterozygousProp+plotDownstreamTotalHeterozygousProp+vertLines+ theme_classic())
# 
# yaxis <- scale_y_continuous(name="Precent of SNPs")
# plotTitle <- ggtitle(paste(outputPrefix, "Saaz vs Frohberg region precent of SNPs", sep=" "))
# plot(plot+plotTitle+xaxis+yaxis+line+plotUpstreamSNPhomozygousProp+plotRegionSNPhomozygousProp+plotDownstreamSNPhomozygousProp+plotUpstreamSNPheterozygousProp+plotRegionSNPheterozygousProp+plotDownstreamSNPheterozygousProp+vertLines+vertLines+ theme_classic())
# 
# 
# 
# 
