library(PopGenome)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

args <- commandArgs(TRUE)
inputName <- args[1] #The folder name where all fasta are stored
chromosomeInfo <- args[2] #"pairwiseDivFile_chromosome.txt"
strainInfo_csv <- args[3]
outputName <- args[4]

chromosome_properties = read.csv(chromosomeInfo,header=TRUE,sep='\t')
chrList <- c(as.character(chromosome_properties$Chr))
chrLength <- c(as.numeric(chromosome_properties$Position_Genome_Start))

list_alignment_files <- list.files(inputName, pattern= 'fasta')
check_strains <- list_alignment_files[1]
listStrains <- c()
for(line in readLines(paste0(inputName,check_strains))){
	if(grepl(">",line)){
		StrainName <- unlist(strsplit(line,">"))[2]
		listStrains <- c(listStrains,StrainName)
	}
}

populations <- list()
counter_Pop <- 1
counterA <- 1
listComparisons <- c()
while(counterA < length(listStrains)){
	counterB <- counterA + 1
	while(counterB <= length(listStrains)){
		assign(paste0(listStrains[counterA],'.VS.',listStrains[counterB]),c(listStrains[counterA],listStrains[counterB]))
		populations[counter_Pop] <-  list(get(paste0(listStrains[counterA],'.VS.',listStrains[counterB]))) #It took me a while to figure this step out
		listComparisons <- c(listComparisons, paste0(listStrains[counterA],'.VS.',listStrains[counterB]))
		counterB <- counterB + 1
		counter_Pop <- counter_Pop + 1
	}
	counterA <- counterA + 1
}

genome <- readData(inputName, include.unknown=FALSE)
genome <- set.populations(genome, populations)
genome <- neutrality.stats(genome)
genome <- diversity.stats(genome, pi=TRUE)
nSites <- genome@n.sites
validSites <- genome@n.valid.sites
piValues <- genome@Pi
regionNames <- genome@region.names

InfoFileName <- paste(paste0("./",outputName), "InfoWindows.txt", sep="_")
con <- file(InfoFileName, open = "wt")
sink(con)
cat("Information regarding data\n")
cat(paste("Median nSites:",median(nSites),"+-",sd(nSites),"\n",sep=" "))
cat(paste("Median validSites:",median(validSites),"+-",sd(validSites),"\n",sep=" "))
sink()
close(con)

output <- data.frame()
for (i in 1:length(regionNames)){
	regionName <- unlist(strsplit(regionNames[i], "[.]"))[1]
	#print(regionName)
	regionSplit <- strsplit(regionName, "_")
	regionChr <- regionSplit[[1]][1]
	regionRange <- strsplit(regionSplit[[1]][2], "-")
	regionStart <- regionRange[[1]][1]
	regionEnd <- regionRange[[1]][2]
	regionAllSites <- nSites[[i]]
	regionValidSites <- validSites[[i]]
	chrIndex <- which(chrList ==regionChr)
	Genome_Pos <- as.numeric(regionStart)+chrLength[chrIndex]
	regionOut <- data.frame(Genome_Pos)
	regionOut$chromosome <- regionChr
	regionOut$Chr_Position <- regionStart
	regionOut$nSites <- regionAllSites
	regionOut$validSites <- regionValidSites
	counterA <- 1
	while(counterA <= length(listComparisons)){
		regionOut[1,4+counterA] <- piValues[i,counterA]
		colnames(regionOut)[4+counterA] <- paste0("rawCounts.",listComparisons[counterA])
		counterA <- counterA + 1
	}
	counterB <- length(colnames(regionOut)) + 1
	counterA <- 1
	while(counterA <= length(listComparisons)){
		regionOut[1,counterB] <- 100*(piValues[i,counterA]/regionValidSites)
		colnames(regionOut)[counterB] <- listComparisons[counterA]
		counterA <- counterA + 1
		counterB <- counterB + 1
	}
	output <- rbind(output, regionOut) 
}
genomeFileName <- paste(paste0("./",outputName), "pairwisePiComp.txt", sep="_")
write.table(output, file=genomeFileName, row.names=FALSE, sep = "\t")

InfoStrains <- read.csv(strainInfo_csv, header = TRUE)
pops2notConsider <- c("Admixture","")
InfoStrains <- subset(InfoStrains, Strains %in% listStrains)
pop_flags <- unique(as.character(subset(InfoStrains, !(Pop %in% pops2notConsider))$Pop))
InfoStrains_popsdf <- subset(InfoStrains, !(Pop %in% pops2notConsider))
print(InfoStrains)

for(i in pop_flags){
  pops2include <- unique(as.character(InfoStrains_popsdf$Pop[InfoStrains_popsdf$Pop != i]))
  for(j in pops2include){
    if(!(exists(paste0("Strains2Compare_",gsub(" ","_",j))))){
      strains2include <- as.character(InfoStrains_popsdf$Strains[InfoStrains_popsdf$Pop == j])
      assign(paste0("Strains2Compare_",gsub(" ","_",j)),strains2include[strains2include != ""])
    }
  }
}

counterA <- length(colnames(output)) + 1
StrainsIndf <- unique(as.character(InfoStrains$Strains))
list_total_comparisons <- c()
listLog2Strains <- c()
for(i in StrainsIndf){
	popStrain <- as.character(InfoStrains$Pop[InfoStrains$Strains == i])
	if(popStrain == "Admixture"){
		listLog2Strains <- c(listLog2Strains,i)
	}
	for(j in pop_flags[pop_flags != popStrain]){
		Strains2Compare_temp <- get(paste0("Strains2Compare_",gsub(" ","_",j)))
		listComparisons_temp <- c()
		for(k in Strains2Compare_temp){
			for(w in listComparisons){
				if(grepl(i,w) & grepl(k,w)){
					listComparisons_temp <- c(listComparisons_temp,w)
				}
			}
		}
		output[,counterA] <- rowMeans(subset(output, select = listComparisons_temp), na.rm=TRUE)
		colnames(output)[counterA] <- paste0(i,".VS.",gsub(" ","_",j),".avg")
		counterA <- counterA + 1
		output[,counterA] <- apply(subset(output, select = listComparisons_temp), 1, sd)
		colnames(output)[counterA] <- paste0(i,".VS.",gsub(" ","_",j),".sd")
		counterA <- counterA + 1
		output[,counterA] <- apply(subset(output, select = listComparisons_temp), 1, min)
		colnames(output)[counterA] <- paste0(i,".VS.",gsub(" ","_",j),".min")
		counterA <- counterA + 1
		list_total_comparisons <- c(list_total_comparisons,paste0(i,".VS.",gsub(" ","_",j)))
	}
}

if(length(listLog2Strains) > 0){
	list_total_comparisons_log2 <- c()
	counterA <- length(colnames(output)) + 1
	for(i in listLog2Strains){
		counterB <- 1
		while(counterB < length(pop_flags)){
			counterC <- counterB + 1
			while(counterC <= length(pop_flags)){
				colnameB <- list_total_comparisons[grepl(i,list_total_comparisons) & grepl(paste0(gsub(" ","_",pop_flags[counterB])),list_total_comparisons)]
				colnameB <- paste0(gsub(" ","_",colnameB),".min")
				colnameC <- list_total_comparisons[grepl(i,list_total_comparisons) & grepl(paste0(gsub(" ","_",pop_flags[counterC])),list_total_comparisons)]
				colnameC <- paste0(gsub(" ","_",colnameC),".min")
				output[,counterA] <- log2(output[[colnameB]]/output[[colnameC]])
				colname_df <- paste0("log.",i,".VS.",gsub(" ","_",pop_flags[counterB]),"_",gsub(" ","_",pop_flags[counterC]))
				list_total_comparisons_log2 <- c(list_total_comparisons_log2,colname_df)
				colnames(output)[counterA] <- colname_df
				counterA <- counterA + 1
				counterB <- counterB + 1
				counterC <- counterC + 1
			}
		}
	}
}

is.na(output) <- sapply(output, is.infinite)
output <- output[with(output, order(Genome_Pos)),]
write.table(output, file=genomeFileName, row.names=FALSE, sep = "\t")


chrs <- unlist(list(unique(output$chromosome)))
chrBreaks <- c()
chrLabel <- c()
lineBreaks <- c()
for (i in 1:length(chrs)) {
	chr_name = toString(chrs[i])
	chr_popG <- subset(output, chromosome==chr_name)
	chrLabel <- append(chrLabel, chr_name)
	lineBreaks <- append(lineBreaks, chr_popG[1,1])
	subLen <- length(chr_popG[,1])
	midpoint <- chr_popG[1,1]+(as.numeric(chr_popG[subLen,3])/2)
	chrBreaks <- append(chrBreaks, midpoint)
}
endPos <- output[length(output[,1]),1]
vertLines <- geom_vline(xintercept = lineBreaks)
line <- geom_abline(intercept=0, slope=0)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

pdf(paste(genomeFileName, "_graphs.pdf", sep="_"), width=14)
#Plots one strain vs one strain
for(i in listStrains){
  comparisons_df <- data.frame(Genome_Pos=numeric(),
                 chromosome=factor(), 
                 Pi_values=numeric(), 
                 Comparison=factor(),
                 stringsAsFactors=FALSE)
  data_strain_temp <- c()
  for(j in listComparisons){
    if(grepl(i,j) & !grepl("rawCounts",j) & !grepl("avg",j) & !grepl("sd",j) & !grepl("min",j)){
      data_strain_temp <- c(data_strain_temp,j)
    }
  }
  for(k in data_strain_temp){
    comparisons_df_temp <- cbind(output[,c(1,2)],output[[k]])
    comparisons_df_temp <- mutate(comparisons_df_temp, Comparison = k)
    colnames(comparisons_df_temp) <- c("Genome_Pos","chromosome","Pi_values","Comparison")
    comparisons_df <- rbind(comparisons_df,comparisons_df_temp)
  }
  comparisons_df$Pi_values[is.na(comparisons_df$Pi_values)] <- -0.1
  plotTitle <- ggtitle(paste0(i," pairwise diversity vs individual Strains"))
  yaxis_breaks <- seq(0,max(comparisons_df$Pi_values), by = 0.5)
  colourCount = length(unique(comparisons_df$Comparison))
  print(ggplot(comparisons_df) +
  geom_point(aes(x=Genome_Pos, y= Pi_values, fill = Comparison, color = Comparison), alpha=0.9) +
  plotTitle +
  theme_classic() +
  scale_x_continuous(breaks=chrBreaks, labels=chrLabel, name="Genome Position", limits=c(0,endPos), expand = c(0,0)) + 
  scale_y_continuous(name="Pairwise diversity (%)", breaks = c(yaxis_breaks,0), limits = c(-0.1,1.5)) +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_color_manual(values = getPalette(colourCount)) +
  line +
  vertLines)
}

#Plots one strain vs populations
for(i in listStrains){
  comparisons_df <- data.frame(Genome_Pos=numeric(),
                 chromosome=factor(), 
                 Average=numeric(),
                 SD=numeric(),
                 ymin_temp=numeric(),
                 ymax_temp=numeric(),
                 Comparison=factor(),
                 stringsAsFactors=FALSE)
  data_strain_temp <- c()
  for(j in pop_flags){
    for(w in list_total_comparisons){
      if(grepl(i,w) & grepl(gsub(" ","_",j),w)){
        data_strain_temp <- c(data_strain_temp,w)
      }
    }
  }
  for(k in data_strain_temp){
    comparisons_df_temp <- cbind(output[,c(1,2)],output[[paste0(k,".avg")]],output[[paste0(k,".sd")]])
    comparisons_df_temp <- mutate(comparisons_df_temp, ymin_temp = comparisons_df_temp[,3]-comparisons_df_temp[,4], ymax_temp = comparisons_df_temp[,3]+comparisons_df_temp[,4], Comparison = k)
    colnames(comparisons_df_temp) <- c("Genome_Pos","chromosome","Average","SD","ymin_temp","ymax_temp","Comparison")
    comparisons_df <- rbind(comparisons_df,comparisons_df_temp)
  }
  comparisons_df$Average[is.na(comparisons_df$Average)] <- -0.1
  comparisons_df$SD[is.na(comparisons_df$SD)] <- 0
  plotTitle <- ggtitle(paste0(i,"Average pairwise diversity"))
  yaxis_breaks <- seq(0,max(comparisons_df$Average), by = 0.5)
  print(ggplot(comparisons_df) +
      geom_point(aes(x=Genome_Pos, y= Average, fill = Comparison, color = Comparison), alpha=0.5) +
      geom_ribbon(aes(x=Genome_Pos, ymin= ymin_temp, ymax = ymax_temp, fill = Comparison, color = Comparison), alpha = 0.3) +
      plotTitle +
      theme_classic() +
      scale_x_continuous(breaks=chrBreaks, labels=chrLabel, name="Genome Position", limits=c(0,endPos), expand = c(0,0)) + 
      scale_y_continuous(name="Pairwise diversity (%)", breaks = c(yaxis_breaks,0)) +
      line +
      vertLines)
}

#Plots log2 admixture strains
for(i in list_total_comparisons_log2){
  comparisons_log_df <- data.frame(Genome_Pos=numeric(),
                 chromosome=factor(), 
                 LogValues=numeric(),
                 stringsAsFactors=FALSE)
  comparisons_log_df_temp <- cbind(output[,c(1,2)],output[[i]])
  comparisons_log_df <- rbind(comparisons_log_df,comparisons_log_df_temp)
  colnames(comparisons_log_df) <- c("Genome_Pos","chromosome","LogValues")
  print(ggplot(comparisons_log_df)+
    geom_line(aes(x=Genome_Pos,y=LogValues, color=i), alpha=1)+
    ggtitle(i) +
      theme_classic() +
      scale_x_continuous(breaks=chrBreaks, labels=chrLabel, name="Genome Position", limits=c(0,endPos), expand = c(0,0)) + 
      line +
      vertLines)
}

print("Done!")