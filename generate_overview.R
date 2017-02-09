#!/usr/bin/env Rscript
library("gplots")
library("graphics")
library("RColorBrewer")
library("MASS")
library("colorRamps")

## arguments from commandline
myarg <- commandArgs(trailingOnly = FALSE)
scriptPath <- grep("--file=", myarg, value=TRUE)
scriptPath <- gsub("--file=", "", scriptPath, fixed=TRUE)
scriptDir <- dirname(scriptPath)
source(file.path(scriptDir, "MethAnR.R"))
argPos <- grep("--args", myarg, fixed = TRUE)
rDir <- as.character(myarg[argPos+1])
outPrefix <- as.character(myarg[argPos+2])
allPositionsFile <- as.character(myarg[argPos+3])
selectedPositionsFile <- as.character(myarg[argPos+4])
gf2dmcFile <- as.character(myarg[argPos+5])
g2gfFile <- as.character(myarg[argPos+6])
g2numFile <- as.character(myarg[argPos+7])
TEGtoTEfile <- as.character(myarg[argPos+8]) # we only care about the second column here

f.print.message <- function(x) { cat("=== generate_overview.R ===", format(Sys.time(), "%Y %b %d %X"), paste0("=== ", x,"\n")) }

######################################################################################
## data reading
######################################################################################

f.print.message("reading data")

myDataFullTable <- read.table(allPositionsFile, sep='\t', header=F,
                              stringsAsFactors=F, row.names=NULL, col.names=c("chrom", "pos"))

myData <- read.table(selectedPositionsFile, sep='\t', header=F,
                     stringsAsFactors=F, row.names=NULL)
colnames(myData) <- paste0("generic_", 1:ncol(myData))
colnames(myData)[1:7] <- c("chrom", "strand", "pos", "meth", "unmeth", "context", "mapping")

rownames(myData) <- paste0(myData$chrom, '_', myData$pos)


######################################################################################
## test different region parameters
######################################################################################

f.print.message("testing region parameters")
res <- f.analyze.regions(myData, rDir, outPrefix, c(seq(5, 20, by = 5), seq(30, 100, by = 10)),
                         c(10,20,30,40,50,75,100,200,300), useLog = TRUE)


######################################################################################
## create regions with a given set of parameters
######################################################################################

#f.print.message("creating regions with minNum=10 and maxDist=25")
#myRegs <- f.create.regions(myData, minNum = 10, maxDist = 25)

# if you want to extract all the cytosines which are within the regions
#myCsInRegs <- f.extract.Cs.within.regions(myData, myRegs)
#myDataSub <- myData[myCsInRegs,]


######################################################################################
## draw a histogram with the distances between significant and random Cs
######################################################################################

f.print.message("drawing distances between Cs")
f.draw.distance.between.DMCs.combined(myData, myDataFullTable, rDir, 4, outPrefix)


######################################################################################
## plot number of cytosines per gene
######################################################################################

f.print.message("loading genes2xyz")
genes2xyzWithTEG <- f.get.genes2xyz.via.files(gf2dmcFile, g2gfFile, g2numFile)
genes2xyzWithoutTEG <- f.remove.tegs.from.genes2xyz(genes2xyzWithTEG, TEGtoTEfile)

f.print.message("plotting number of Cs per gene")
f.plot.number.DMCs.per.gene(genes2xyzWithTEG, rDir, outPrefix)
f.plot.number.DMCs.per.gene(genes2xyzWithoutTEG, rDir, paste0(outPrefix, "_withoutTEOR"))

f.print.message("FINISHED")
