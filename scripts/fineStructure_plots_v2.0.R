source("~/software/scripts/sources/FineStructureLibrary.R")
library(optparse)
library(ape)
library(paran)
library(XML)
library(psych)

 
option_list = list(
make_option(c("-p", "--prefix"), action="store", type="character", default=NULL, 
  help="Name of the output name in finestructure, for example Smik_AllGenomes_fs_linked_hap.chunkcounts.out, your names is Smik_AllGenomes_fs",metavar="character"),
make_option(c("-s", "--state"), action="store", type="character", default=NULL, 
  help="output file name linked/unlinked; if you used the recombination map it will result in linked results",metavar="character"),
make_option(c("-I", "--strainInfo_csv"), action="store", type="character", default=NULL, 
  help="ADDED by me to take strainINFO information",metavar="character"),
make_option(c("-i", "--inputfolder"), action="store", type="character", default=NULL, 
  help="Where the fineStructure outputs are stored",metavar="character"),
make_option(c("-N", "--prefix_NameInFasta"), action="store", type="character", default="", 
  help="If your strain names contain an addition, we can try to remove them,  [default= %default].",metavar="character"),
make_option(c("-a", "--admixedList"), action="store", type="character", default=NULL, 
  help="If you have a list of known admixture strains, but the functionallity of this is still not implemented",metavar="character")
)
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

prefix <- opt$prefix
state <- opt$state
strainInfo_csv <- opt$strainInfo_csv
inputfolder <- opt$inputfolder
prefix_NameInFasta <- opt$prefix_NameInFasta
#admixedList <- opt$admixedList

strainInfo_csv <- read.csv(strainInfo_csv,header = TRUE)

## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

### Define our input files
chunkfile<-paste(paste0(inputfolder,prefix), state , "hap.chunkcounts.out", sep="_") ## chromopainter chunkcounts file
mcmcfile<-paste(paste0(inputfolder,prefix), state, "hap_mcmc_run0.xml", sep="_") ## finestructure mcmc file
mcmcfileR1<-paste(paste0(inputfolder,prefix), state, "hap_mcmc_run1.xml", sep="_") ## finestructure mcmc file
treefile<-paste(paste0(inputfolder,prefix), state, "hap_tree_run0.xml", sep="_") ## finestructure tree file
treefileR1<-paste(paste0(inputfolder,prefix), state, "hap_tree_run1.xml", sep="_") ## finestructure tree file

list_text2modify <- c(chunkfile,mcmcfile,mcmcfileR1,treefile,treefileR1)
#Added in case we want to remove the prefix in the fasta IDs
if(prefix_NameInFasta != ""){
	for(i in list_text2modify){
	tx_temp <- readLines(i)
	tx2_temp <- gsub(pattern=prefix_NameInFasta, replace ="", x=tx_temp)
	writeLines(tx2_temp, con=i)
	}
}

## Additional files that you can extract from finestructure
mappopchunkfile<-paste(paste0(inputfolder,prefix), state, "hap_mapstate.csv", sep="_") # population-by-population chunkcount file for the populations used in the MAP (i.e tree)
mappopchunkfileR1<-paste(paste0(inputfolder,prefix), state, "hap_mapstateR1.csv", sep="_")
print(mappopchunkfile)

system( paste("fs fs -X -Y -e X2",chunkfile,treefile,mappopchunkfile) )
system( paste("fs fs -X -Y -e X2",chunkfile,treefileR1,mappopchunkfileR1) )
meancoincidencefile<-paste(paste0(inputfolder,prefix), state, "hap_meancoincidence.csv", sep="_")
meancoincidencefileR1<-paste(paste0(inputfolder,prefix), state, "hap_meancoincidenceR1.csv", sep="_") # pairwise coincidence, .i.e. proportion of MCMC files where individuals are found in the same 
system( paste("fs fs -X -Y -e meancoincidence",chunkfile,mcmcfile,meancoincidencefile) )
system( paste("fs fs -X -Y -e meancoincidence",chunkfile,mcmcfileR1,meancoincidencefileR1) )
## there are ways of generating these within R but are either slower or more annoying - its your call how you do it

###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 

###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile)
mcmcxmlR1<-xmlTreeParse(mcmcfileR1) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml)
mcmcdataR1<-as.data.frame.myres(mcmcxmlR1) ## convert this into a data frame

###### READ IN THE TREE FILES
treexml<-xmlTreeParse(treefile)
treexmlR1<-xmlTreeParse(treefileR1) ## read the tree as xml format
ttree<-extractTree(treexml)
ttreeR1<-extractTree(treexmlR1) ## extract the tree into ape's phylo format
## If you dont want to plot internal node labels (i.e. MCMC posterior assignment probabilities)
## now is a good time to remove them via:
#     ttree$node.label<-NULL
## Will will instead remove "perfect" node labels
ttree$node.label[ttree$node.label=="1"] <-""
ttreeR1$node.label[ttreeR1$node.label=="1"] <-""
## And reduce the amount of significant digits printed:
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
ttreeR1$node.label[ttreeR1$node.label!=""] <-format(as.numeric(ttreeR1$node.label[ttreeR1$node.label!=""]),digits=2)

tdend<-myapetodend(ttree,factor=1) # convert to dendrogram format
tdendR1<-myapetodend(ttreeR1,factor=1)

####################################
## PLOT 1: RAW DENDROGRAM PLOT
pdf(file=paste(paste0(inputfolder,prefix), state, "FullDendrogram.pdf", sep="_"),height=6,width=14)
par(mar=c(6,0,2,0),mfrow=c(1,1))
fs.plot.dendrogram(tdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=1),edgePar=list(p.lwd=0,t.srt=90,t.off=-0.5),axes=F)
dev.off()

pdf(file=paste(paste0(inputfolder,prefix), state, "FullDendrogramR1.pdf", sep="_"),height=6,width=14)
par(mar=c(6,0,2,0),mfrow=c(1,1))
fs.plot.dendrogram(tdendR1,horiz=FALSE,nodePar=list(cex=0,lab.cex=1),edgePar=list(p.lwd=0,t.srt=90,t.off=-0.5),axes=F)
dev.off()

## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstateR1<-extractValue(treexmlR1,"Pop")
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
mapstatelistR1<-popAsList(mapstateR1)

popnames<-lapply(mapstatelist,NameSummary)
popnamesR1<-lapply(mapstatelistR1,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
popnamesplotR1<-lapply(mapstatelistR1,NameMoreSummary)
names(popnames)<-popnamesplot # for nicety only
names(popnamesR1)<-popnamesplotR1
names(popnamesplot)<-popnamesplot # for nicety only
names(popnamesplotR1)<-popnamesplotR1

#HERE there is a step added by Quinn to define pops, but first let's check the second chunk

total_pops <- unique(as.character(strainInfo_csv$Pop))
pops<- c()
for(i in total_pops){
  assign(gsub(" ","_",i),as.character(strainInfo_csv$Strains[strainInfo_csv$Pop == i]))
  pops[[gsub(" ","_",i)]] <- get(gsub(" ","_",i))
}

###
#popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
#popdend<-makemydend(tdend,mapstatelist)
popdend<-makemydend(tdend,pops)
popdendR1<-makemydend(tdendR1,pops)
#popdend<-fixMidpointsComplete(popdend) # needed for obscure dendrogram reasons
popdend<-fixMidpointMembers(popdend) #This does exist but I don't know what it does
popdendR1<-fixMidpointMembers(popdendR1)

########################
## PLOT 2: population tree
pdf(file=paste(paste0(inputfolder,prefix), state, "PopulationDendrogram.pdf", sep="_"),height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
fs.plot.dendrogram(popdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=1,las=2),edgePar=list(p.lwd=0,t.srt=0,t.off=0.3),yaxt="n",height=0.5,dLeaf=0.2)
dev.off()

pdf(file=paste(prefix, state, "PopulationDendrogramR1.pdf", sep="_"),height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
fs.plot.dendrogram(popdendR1,horiz=FALSE,nodePar=list(cex=0,lab.cex=1,las=2),edgePar=list(p.lwd=0,t.srt=0,t.off=0.3),yaxt="n",height=0.5,dLeaf=0.2)
dev.off()

########################
## PAIRWISE COINCIDENCES
fullorder<-labels(tdend)
fullorderR1<-labels(tdendR1) # the order according to the tree
mcmcmatrixraw<-as.matrix(read.csv(meancoincidencefile,row.names=1))
mcmcmatrixrawR1<-as.matrix(read.csv(meancoincidencefileR1,row.names=1)) # read in the pairwise coincidence file we created earlier
mcmcmatrixR1<-mcmcmatrixrawR1[fullorderR1,fullorderR1] 
mcmcmatrix<-mcmcmatrixraw[fullorder,fullorder] 
mapstatematrix<-groupingAsMatrix(pops)[fullorder,fullorder]
mapstatematrixR1<-groupingAsMatrix(pops)[fullorderR1,fullorderR1] # map state for reference

#########################
## PLOT 3: Pairwise coincidence, showing the MAP state

#source("FinestructureLibrary.R")
pdf(file=paste(paste0(inputfolder,prefix), state, "PairwiseCoincidence.pdf", sep="_"),height=12,width=12)
plotFinestructure(mcmcmatrix,dimnames(mcmcmatrix)[[1]],dend=tdend,optpts=mapstatematrix,cex.axis=1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

pdf(file=paste(paste0(inputfolder,prefix), state, "PairwiseCoincidenceR1.pdf", sep="_"),height=12,width=12)
plotFinestructure(mcmcmatrixR1,dimnames(mcmcmatrixR1)[[1]],dend=tdendR1,optpts=mapstatematrixR1,cex.axis=1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1))
dev.off()

########################
## COANCESTRY MATRIX

datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix
datamatrixR1<-dataraw[fullorderR1,fullorderR1]

tmatmax<-max(datamatrix)+sd(datamatrix) # cap the heatmap
tmpmat<-datamatrix 
tmpmatR1<-datamatrixR1 
tmpmat[tmpmat>tmatmax]<-tmatmax # 
tmpmatR1[tmpmatR1>tmatmax]<-tmatmax
pdf(file=paste(paste0(inputfolder,prefix), state, "Coancestry.pdf", sep="_"),height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

pdf(file=paste(paste0(inputfolder,prefix), state, "CoancestryR1.pdf", sep="_"),height=12,width=12)
plotFinestructure(tmpmatR1,dimnames(tmpmatR1)[[1]],dend=tdendR1,cols=some.colorsEnd,cex.axis=1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

## Population averages
popmeanmatrix<-getPopMeanMatrix(datamatrix,pops)
popmeanmatrixR1<-getPopMeanMatrix(datamatrixR1,pops)

tmatmax<-max(datamatrix)+sd(datamatrix) # cap the heatmap (we need to increase)
tmpmat<-popmeanmatrix
tmpmatR1<-popmeanmatrixR1
tmpmat[tmpmat>tmatmax]<-tmatmax
tmpmatR1[tmpmatR1>tmatmax]<-tmatmax# 
pdf(file=paste(paste0(inputfolder,prefix), state, "PopAveragedCoancestry.pdf", sep="_"),height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

pdf(file=paste(paste0(inputfolder,prefix), state, "PopAveragedCoancestryR1.pdf", sep="_"),height=12,width=12)
plotFinestructure(tmpmatR1,dimnames(tmpmat)[[1]],dend=tdendR1,cols=some.colorsEnd,cex.axis=1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()

### Useful tricks with labels
#mappopcorrectorder<-NameExpand(labels(popdend))
#mappopcorrectorderR1<-NameExpand(labels(popdendR1))
#mappopsizes<-sapply(mappopcorrectorder,length)
#mappopsizesR1<-sapply(mappopcorrectorderR1,length)
#labellocs<-PopCenters(mappopsizes)
#labellocsR1<-PopCenters(mappopsizesR1)
#labelcols<-c(2,2,3,3,4,4,1,1,1,5,6,6,6,7,8,8,2) # different label colours allow clearer identification of individuals, too
#labelcrt=20

#pdf(file=paste(paste0(inputfolder,prefix), state, "Coancestry2.pdf", sep="_"),height=12,width=12)
#plotFinestructure(tmpmat,labelsx=labels(popdend),labelsatx=labellocs,crt=labelcrt,dend=tdend,text.col=labelcols,cex.axis=1.0,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
#dev.off()

#pdf(file=paste(paste0(inputfolder,prefix), state, "Coancestry2R1.pdf", sep="_"),height=12,width=12)
#plotFinestructure(tmpmatR1,labelsx=labels(popdendR1),labelsatx=labellocsR1,crt=labelcrt,dend=tdendR1,text.col=labelcols,cex.axis=1.0,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
#dev.off()

####################################
## PCA Principal Components Analysis
pcares<-mypca(dataraw)
eigenvalues <- pcares$values
# For figuring out how many PCs are important; see Lawson & Falush 2012
# You need packages GPArotation and paran
tmap<-optimalMap(dataraw)
thorn<-optimalHorn(dataraw)
c(tmap,thorn) # 11 and 5. Horn typically underestimates, Map is usually better
pcapops<-getPopIndices(rownames(dataraw),pops)
pcanames<-rownames(dataraw)
rcols<-rainbow(max(pcapops))
colors <- c()
for(i in pcapops){
  colors <- c(colors,rcols[i])
}

pdf(paste(paste0(inputfolder,prefix), state, "PCA1-3.pdf",sep="_"), height=6,width=12)
par(mfrow=c(1,2))
i=1
for(j in (i+1):3) {
  plot(pcares$vectors[,i],pcares$vectors[,j],col="black",bg=colors,xlab=paste("PC",i,as.integer(eigenvalues[i]/sum(eigenvalues)*100),"%"),ylab=paste("PC",j,as.integer(eigenvalues[j]/sum(eigenvalues)*100),"%"),main=paste("PC",i,"vs",j),pch=21,cex=1)
}
par(mfrow=c(4,3))
for(i in 1:4) for(j in (i+1):5) {
  plot(pcares$vectors[,i],pcares$vectors[,j],col=rcols[pcapops],xlab=paste("PC",i,as.integer(eigenvalues[i]/sum(eigenvalues)*100),"%"),ylab=paste("PC",j,as.integer(eigenvalues[j]/sum(eigenvalues)*100),"%"),main=paste("PC",i,"vs",j),pch=16)
  text(pcares$vectors[,i],pcares$vectors[,j],labels=pcanames,col=rcols[pcapops],cex=1,pos=1)
}
dev.off()

#########################
## CHROMOPAINTER Section is turne OFF add in the future see chooseFineStrPlots_modifiedPCA.R for future additions
