#!/usr/bin/env Rscript
library("plotrix")
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
sigStatsFile <- as.character(myarg[argPos+3])
species <- as.character(myarg[argPos+4])

f.print.message <- function(x) { cat("=== generate_mapping_stats_figures.R ===", format(Sys.time(), "%Y %b %d %X"), paste0("=== ", x,"\n")) }

######################################################################################
## check mapping stats
######################################################################################

f.print.message("loading mapping stats")
mappingStats <- read.table(sigStatsFile, sep='\t', header=T, row.names=NULL, fill=T, stringsAsFactors = FALSE)
mappingStats <- mappingStats[,c("context", "feature", "methylatedCs", "unmethylatedCs", "methylatedCoverage", "unmethylatedCoverage")]
colnames(mappingStats) <- c("context", "feature", "meCnum", "deCnum", "methCov", "demethCov")
mappingStats$totalC <- mappingStats$meCnum + mappingStats$deCnum
mappingStats$perc <- round(mappingStats$totalC/sum(mappingStats$totalC)*100,3)

myContexts <- intersect(c("CG", "CpG"), unique(mappingStats[,"context"]))
if (length(myContexts) == 0) {
  myContexts <- c("CG", "CHG", "CHH")
} else {
  myContexts <- c(myContexts, "CHG", "CHH")
}


if (species == "At") {
  baseLine <- c(CG=13,CpG=13,CHG=14,CHH=73) # Arabidopsis thaliana
} else if (species == "Mp") {
  baseLine <- c(CG=18,CpG=18,CHG=17,CHH=65) # Marchantia polymorpha
} else {
  f.print.message("WARNING: THE SPECIES YOU SPECIFIED DOES NOT HAVE A DEFAULT")
  baseLine <- c(CG=0,CpG=0,CHG=0,CHH=0) # THATS ORGANISM-SPECIFIC
}

f.print.message("plotting pie charts")
svg(file.path(rDir, paste0(outPrefix, "_dmcPieCharts.svg")), height = 10, width = 10)
layout(matrix(1:4, ncol = 2))
sry <- aggregate(cbind(totalC, meCnum, deCnum) ~ context, data = mappingStats, sum)
if (sum(baseLine)>0) { pie(baseLine[sry$context], labels = sry$context, col = c("gray10", "gray40", "gray80"), main = "Cs") }
pie(sry$totalC, labels = sry$context, col = c("gray10", "gray40", "gray80"), main = "all Cs")
if (sum(sry$meCnum) > 0) { pie(sry$meCnum, labels = sry$context, col = c("gray10", "gray40", "gray80"), main = "methylated Cs") }
if (sum(sry$deCnum) > 0) { pie(sry$deCnum, labels = sry$context, col = c("gray10", "gray40", "gray80"), main = "unmethylated Cs") }
dev.off()
write.csv(sry, file.path(rDir, paste0(outPrefix, "_pieChartNumbers.csv")), quote = FALSE)

f.print.message("preparing for spider graphs")
contextColors <- data.frame(
  CG = brewer.pal(3,"YlOrBr"),
  CHG = brewer.pal(3,"PuBuGn"),
  CHH = brewer.pal(3,"BuGn"),
  row.names = c("totalC", "meCnum", "deCnum"), stringsAsFactors = FALSE
); names(contextColors) <- myContexts
segPlots <- paste(rep(myContexts, each = 2), rep(c("meCnum", "deCnum"), 3), sep = '_')
segCols <- unique(mappingStats$feature)
numPlots <- matrix(0, nrow = length(segPlots), ncol = length(segCols), dimnames = list(segPlots, segCols))
for (ctxt in myContexts) {
  temp <- subset(mappingStats, context == ctxt)
  for (pc in c("meCnum", "deCnum")) {
    numPlots[paste(ctxt, pc, sep = '_'), temp$feature] <- temp[[pc]]
  }
}
numPlots <- as.data.frame(numPlots)

# numPlots$transposon <- numPlots$transposable_element + numPlots$transposable_element_gene
# # rearrange manually for nice order and plot spider graphs (see colnames(numPlots))
# numPlots <- numPlots[,c("intergenic", "transposon", "five_flank", "five_prime_UTR", "exon",
#                         "splice", "three_prime_UTR", "three_flank", "pseudogene",
#                         "miRNA_primary_transcript")]
# labs <- c("Intergenic", "Transposon", "5'Upstream", "5'UTR", "Exon", "Intron",
#           "3'UTR", "3'Downstream", "Pseudogene", "miRNA")

f.print.message("summarizing annotation")
if (species == "At") {
  transposonOrRepeatTypes <- list(
    transposon = c("transposable_element", "transposable_element_gene"),
    pseudogene = c("pseudogene"), # yes - that's not a repeat or transposon :)
    miRNAgene = c("miRNA_primary_transcript") # this as well not
  )
} else if (species == "Mp") {
  # ignore: rRNA, snRNA
  colnames(numPlots) <- gsub("^transposable_element_or_repeat_", "", colnames(numPlots))
  transposonOrRepeatTypes <- list(
    DNAtransp = grep("^DNA", colnames(numPlots), value = TRUE),
    LTRtransp = grep("^LTR", colnames(numPlots), value = TRUE),
    LINEtransp = grep("^LINE", colnames(numPlots), value = TRUE),
    SINEtransp = grep("^SINE", colnames(numPlots), value = TRUE),
    RCtransp = grep("^RC", colnames(numPlots), value = TRUE),
    repeats = c("Satellite", "Satellite_telomeric", "Simple_repeat")#,
    #uknRepeats = c("unClassRep", "Unknown")
  )
} else {
  f.print.message("WARNING: THE SPECIES YOU SPECIFIED DOES NOT HAVE A DEFAULT - using Mp")
  colnames(numPlots) <- gsub("^transposable_element_or_repeat_", "", colnames(numPlots))
  transposonOrRepeatTypes <- list(
    DNAtransp = grep("^DNA", colnames(numPlots), value = TRUE),
    LTRtransp = grep("^LTR", colnames(numPlots), value = TRUE),
    LINEtransp = grep("^LINE", colnames(numPlots), value = TRUE),
    SINEtransp = grep("^SINE", colnames(numPlots), value = TRUE),
    RCtransp = grep("^RC", colnames(numPlots), value = TRUE),
    repeats = c("Satellite", "Satellite_telomeric", "Simple_repeat")#,
    #uknRepeats = c("unClassRep", "Unknown")
  )
}

for (toMerge in names(transposonOrRepeatTypes)) {
  if  (sum(transposonOrRepeatTypes[[toMerge]] %in% colnames(numPlots)) == 0) {
    f.print.message(paste0("INFO: skipping, ", toMerge))
    next
  }
  temp <- numPlots[,transposonOrRepeatTypes[[toMerge]]]
  if (is.vector(temp)) {
    numPlots[[toMerge]] <- temp
  } else {
    numPlots[[toMerge]] <- apply(temp, 1, sum)
  }
}

# rearrange manually for nice order and plot spider graphs (see colnames(numPlots))
TORPs <- names(transposonOrRepeatTypes)
if (length(TORPs) > 1) {
  temp <- floor(length(TORPs)/2)
  firstPart <- TORPs[1:temp]
  secondPart <- TORPs[(temp+1):length(TORPs)]
} else {
  firstPart <- c()
  secondPart <- TORPs
}
niceOrder <- c("Intergenic", firstPart, "5'Upstream", "5'UTR", "Exon", "Intron", "3'UTR", "3'Downstream", secondPart)
names(niceOrder) <- c("intergenic", firstPart, "five_flank", "five_prime_UTR", "exon", "splice", "three_prime_UTR", "three_flank", secondPart)
niceOrder <- niceOrder[intersect(names(niceOrder), colnames(numPlots))]
numPlots <- numPlots[,names(niceOrder)]
labs <- niceOrder

f.print.message("plotting spider graphs")
percPlots <- numPlots/apply(numPlots, 1, sum)*100
acrossPercPlots <- numPlots/sum(numPlots)*100
svg(file.path(rDir, paste0(outPrefix, "_spiderPlots.svg")), height = 15, width = 10)
layout(matrix(1:9, ncol = 3, byrow = TRUE))
for (ctxt in myContexts) {
  cR <- paste(ctxt, c("meCnum", "deCnum"), sep = '_')
  meCol <- contextColors["meCnum", ctxt]
  deCol <- contextColors["deCnum", ctxt]
  percMax <- max(ceiling(percPlots[cR,]/10), na.rm = TRUE)*10
  if (percMax < 10) {percMax <- 10}
  acrossPercMax <- max(ceiling(acrossPercPlots[cR,]/10), na.rm = TRUE)*10
  if (acrossPercMax < 10) {acrossPercMax <- 10}
  radial.plot(percPlots[cR,], labels = labs, main = ctxt, radial.lim = seq(0, percMax,10), 
              start = pi/2, line.col = c(meCol,deCol), poly.col = c(meCol, NA),
              rp.type = "p", show.grid.labels = 3, lwd = 3, point.symbols = 16,
              point.col = c(meCol,deCol), show.centroid = FALSE)
  radial.plot(acrossPercPlots[cR,], labels = labs, main = ctxt, radial.lim = seq(0, acrossPercMax,10),
              start = pi/2, line.col = c(meCol,deCol), poly.col = c(meCol, NA),
              rp.type = "p", show.grid.labels = 3, lwd = 3, point.symbols = 16,
              point.col = c(meCol,deCol), show.centroid = FALSE)
  radial.plot(numPlots[cR,], labels = labs, main = ctxt, 
              start = pi/2, line.col = c(meCol,deCol), poly.col = c(meCol, NA),
              rp.type = "p", show.grid.labels = 3, lwd = 3, point.symbols = 16,
              point.col = c(meCol,deCol), show.centroid = FALSE)
}
dev.off()
write.csv(numPlots, file.path(rDir, paste0(outPrefix, "_spiderGraphNumbers.csv")), quote = FALSE)
f.print.message("FINISHED")
