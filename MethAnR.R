#' TITLE
#'@param
#'@return
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export


###########################################################################################
####### 
###########################################################################################


###########################################################################################
####### general things
###########################################################################################

#' Revert histogram data (vector)
#'@note this function takes a vector of two entries
#'1: the value (for example distance)
#'2: the occurence of the value
#'it will then return a vector that holds the value as many times as it occurs
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.reverse.histogram.data <- function(x) {
  return(rep(x[1], x[2]))
}

#' Revert histogram data (matrix)
#'@note this function takes a vector of two entries
#'1: the value (for example distance)
#'2: the occurence of the value
#'it will then return a vector that holds the value as many times as it occurs
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.reverse.histogram.data.matrix <- function(x) {
  return(unlist(apply(x, 1, f.reverse.histogram.data)))
}

#' Draw a histogram (frequency or density)
#'@param x values
#'@param xName x-axis label
#'@param useFreq plot frequencies instead of densities
#'@param doNotShowSummary omit printing mean/median/sd
#'@param useLog log2-transform the data
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.histogram <- function(x, y = c(), xName = "", useFreq = FALSE, doNotShowSummary = FALSE, useLog = FALSE, ...) 
{
  if (length(x) == length(y)) {
    toPlot <- list(x = x, y = y)
    class(toPlot) <- "density"
  } else {
    toPlot <- density(x)
  }
  if (useFreq) { 
    toPlot$y <- toPlot$y * length(x)
  }
  if (useLog) {
    if (useFreq) { toPlot$y <- log2(toPlot$y+1) }
    else { cat("log only possible for useFreq = TRUE\n") }
  }
  par(mar = c(5,4,1,1))
  plot(
    toPlot,
    #main = "",
    bty = "n",
    xaxs = "r",
    yaxs = "r",
    xlab = xName,
    ylab = "",
    las = 1,
    cex = 0.4,
    tck = 0.01,
    ...
  )
  if (!doNotShowSummary) {
    text(
      par("usr")[1] + par("usr")[2]/15,
      par("usr")[4] - par("usr")[4]/15,
      adj = c(0,1),
      labels = paste(c("median:", "mean:", "sd:"),
                     c(round(median(x), digits = 3),
                       round(mean(x), digits = 3),
                       round(sd(x), digits = 3)),
                     collapse = "\n", sep=" ")
    )
  }
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.genetodensity <- function(x,y) {
  require("MASS")
  est <- kde2d(x, y, n = 50) # needs MASS - calculates the two dimensional density
  if (sum(is.na(est$z)) > 0) {
    est <- kde2d(x, y, n = 50, h = c(bw.nrd0(x), bw.nrd0(y))) # needs MASS - calculates the two dimensional density
  }
  dvec <- as.vector(est$z) # note: as.vector takes columnwise, rows correspond to the value of x, columns to the value of y - the y index counts therefore 50 times
  xind <- findInterval(x, est$x) # find the intervals to which a gene belongs
  yind <- findInterval(y, est$y) # find the intervals to which a gene belongs
  vecind <- (yind-1)*50 + xind # this is a vector with indices for dvec. 
  dens <- dvec[vecind]
  names(dens) <- names(x)
  return(dens)
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.genetodensitycolor <- function(x,y) {
  require("MASS")
  require("colorRamps")
  ## NOTE would not work if two spots have exactly the same density or?
  ## important for a nice picture is that one removes the grid points with no gene (the ones with the lowest density)
  grids <- 100
  out <- list()
  est <- kde2d(x, y, n = grids) # needs MASS - calculates the two dimensional density
  if (sum(is.na(est$z)) > 0) {
    est <- kde2d(x, y, n = grids, h = c(bw.nrd0(x), bw.nrd0(y))) # needs MASS - calculates the two dimensional density
  }
  dvec <- as.vector(est$z) # note: as.vector takes columnwise, rows correspond to the value of x, columns to the value of y - the y index counts therefore grids times
  xind <- findInterval(x, est$x) # find the intervals to which a gene belongs
  yind <- findInterval(y, est$y) # find the intervals to which a gene belongs
  vecind <- (yind-1)*grids + xind # this is a vector with indices for dvec. 
  dens <- dvec[vecind]
  dens_with_genes <- unique(dens)
  names(dens) <- names(x)
  out$dens <- dens
  out$denschar <- as.vector(dens, mode = "character")
  colgrad <- blue2red(length(dens_with_genes))
  #colgrad <- rainbow(length(dens_with_genes), start = 0, end = 1, alpha = 0.8)
  #colgrad <- heat.colors(length(dens_with_genes), alpha = 0.8)
  #colgrad <- terrain.colors(length(dens_with_genes), alpha = 0.8)
  #colgrad <- topo.colors(length(dens_with_genes), alpha = 0.8)
  #colgrad <- cm.colors(length(dens_with_genes), alpha = 0.8)
  sdens_with_genes <- sort(dens_with_genes, decreasing = FALSE)
  names(colgrad) <- sdens_with_genes
  out$cols <- colgrad
  genecols <- colgrad[as.vector(dens, mode = "character")]
  names(genecols) <- names(x)
  out$genecols <- genecols
  return(out)
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.yellowblueblack <- function(x) {
  #rg <- approx(c(0, 0.5, 1), c(1, 1/3, 0), n = x)$y
  #b <- approx(c(0, 0.5, 1), c(2/3, 2/3, 0), n = x)$y
  rg <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 1, 0), n = x)$y
  return(rgb(rg, rg, b))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.yellowblackblue <- function(x) {
  #rg <- approx(c(0, 1, 0.5), c(1, 1/3, 0), n = x)$y
  #b <- approx(c(0, 1, 0.5), c(2/3, 2/3, 0), n = x)$y
  rg <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 0, 1), n = x)$y
  return(rgb(rg, rg, b))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.yellowredblue <- function(x) {
  r <- approx(c(0, 0.5, 1), c(1, 1, 0), n = x)$y
  g <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 0, 1), n = x)$y
  return(rgb(r, g, b))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.yellowblack <- function(x) {
  rg <- approx(c(0, 1), c(1, 0), n = x)$y
  b <- approx(c(0, 1), c(0, 0), n = x)$y
  return(rgb(rg, rg, b))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.yellowredblack <- function(x) {
  r <- approx(c(0, 0.5, 1), c(1, 1, 0), n = x)$y
  g <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 0, 0), n = x)$y
  return(rgb(r, g, b))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.redwhiteblack <- function(x) {
  r <- approx(c(0, 0.5, 1), c(1, 1, 0), n = x)$y
  g <- approx(c(0, 0.5, 1), c(0, 1, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 1, 0), n = x)$y
  return(rgb(r, g, b))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.blackblueyellow <- function(x) {
  return(rev(f.yellowblueblack(x)))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.blueblackyellow <- function(x) {
  return(rev(f.yellowblackblue(x)))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.blueredyellow <- function(x) {
  return(rev(f.yellowredblue(x)))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.blackyellow <- function(x) {
  return(rev(f.yellowblack(x)))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.blackredyellow <- function(x) {
  return(rev(f.yellowredblack(x)))
}

#' undocumented internal function
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.blackwhitered <- function(x) {
  return(rev(f.redwhiteblack(x)))
}

###############################################################################
## plot position/region densities along chromosomes
###############################################################################

#' plot position (e.g. cytosine) density along the (Arabidopsis) genome
#'@param dmcTable a data frame with the positions ("chrom" and "pos" are required)
#'@param rDir the directory in which the results will be stored
#'@param outfile file name without path and extension (default is "posDensityAlongChromosomes")
#'@return NULL
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.plot.pos.density <- function(dmcTable, rDir, outfile = "posDensityAlongChromosomes") {
  #centromeres <- list(Chr1 = 15088987/1e6, Chr2 = 3608426.5/1e6, Chr3 = 13591999.5/1e6, Chr4 = 3956518.5/1e6, Chr5 = 11742754.5/1e6)
  centromeres <- c(15088987/1e6, 3608426.5/1e6, 13591999.5/1e6, 3956518.5/1e6, 11742754.5/1e6)
  cenY <- 5
  ylimit <- c(0,100)
  xlimit <- c(-0.5,30.427671)
  xAchse <- as.character(rep(0:30, each = 2))
  xAchse[seq(2, length(xAchse), 2)] <- ""
  xAchse <- xAchse[1:(length(xAchse)-1)]
  outfileSVG <- paste0(outfile, ".svg")
  outfilePNG <- paste0(outfile, ".png")
  svg(filename = file.path(rDir, outfileSVG), pointsize = 15, width = 23.4, height = 16.6)
  par(oma=c(6,1,1,1))
  layout(mat = c(1:5), width = rep(1,5), height = rep(2,5))
  for (i in 1:5){
    par(mar = c(1,1,1,1))
    subDMC <- subset(dmcTable, chrom == i)
    plot(0, 0, type="n", ylim=ylimit, xlim=xlimit, ylab="", xlab ="",yaxt="n",xaxt="n",cex.lab = 1.5, cex.axis = 1.5,bty="n", pch=1, col = "#00ced11e", lwd = 1)
    #sapply(subDMC$pos/1e6, function(x) lines(c(x,x),c(50,70), col = "#00ced11e", lwd = 1))
    den <- density(subDMC$pos/1e6, bw = 0.01)
    den$y <- den$y*(80/max(den$y))
    lines(den$x, den$y, col = "red")
    # other thins
    points(centromeres[i],cenY,pch=16,cex=3,bg="black")
    mtext(paste0("Chr", i),side=4,cex=1,adj=0.5,padj=-3)
    axis(1, at=seq(0,30,by=0.5),labels=xAchse,col.axis="black",las=0, cex.axis = 1.5, cex.lab = 1.5)
    if (i == 5) { mtext("position (Mb)",side=1,cex=1,adj=0.5,padj=4) }
    #if (i == 2) { legend("right", title = "", legend = c(genotypes,"shared DMC loci", "shared DMCs", "centromere"), col = c(unlist(genotypeCol[genotypes]),"black","black","black"), pch = c(16,16,16,4,16), bty = 'n', cex =1.5, inset = 0.05, pt.cex = 1.5) }
  }
  dev.off()
  system(paste0("rsvg-convert -a -d 300 -p 300 ", file.path(rDir, outfileSVG)," > ", file.path(rDir, outfilePNG)))
  return(NULL)
}

#' plot region density along the (Arabidopsis) genome
#'@param regTable a data frame with the regions ("chrom", "start", "size" and "density" are required)
#'@param rDir the directory in which the results will be stored
#'@param outfile file name without extension (default is "regionDensityAlongChromosomes")
#'@return NULL
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.plot.region.density <- function(regTable, rDir, outfile = "regionDensityAlongChromosomes") {
  #centromeres <- list(Chr1 = 15088987/1e6, Chr2 = 3608426.5/1e6, Chr3 = 13591999.5/1e6, Chr4 = 3956518.5/1e6, Chr5 = 11742754.5/1e6)
  centromeres <- c(15088987/1e6, 3608426.5/1e6, 13591999.5/1e6, 3956518.5/1e6, 11742754.5/1e6)
  cenY <- 5
  ylimit<-c(0,100)
  xlimit<-c(-0.5,30.427671)
  xAchse <- as.character(rep(0:30, each = 2))
  xAchse[seq(2, length(xAchse), 2)] <- ""
  xAchse <- xAchse[1:(length(xAchse)-1)]
  regTable$pos <- regTable$start+regTable$size/2
  outfileSVG <- paste0(outfile, ".svg")
  outfilePNG <- paste0(outfile, ".png")
  svg(filename = file.path(rDir, outfileSVG), pointsize = 15, width = 23.4, height = 16.6)
  par(oma=c(6,1,1,1))
  layout(mat = c(1:5), width = rep(1,5), height = rep(2,5))
  for (i in 1:5){
    par(mar = c(1,1,1,1))
    subTable <- subset(regTable, chrom == i)
    colGrad <- f.blackredyellow(15) # 0 to 1.5 in 0.1
    colBins <- seq(0, 1.5, by = 0.1)
    plot(0, 0, type="n", ylim=ylimit, xlim=xlimit, ylab="", xlab ="",yaxt="n",xaxt="n",cex.lab = 1.5, cex.axis = 1.5,bty="n", pch=1, col = "#00ced10f", lwd = 1)
    sapply(subTable$pos/1e6, function(x) lines(c(x,x),c(10,90), col = "#00ced164", lwd = 1))
    #xVals <- subTable$pos/1e6
    #yVals <- rep(85, nrow(subTable))
    #pCols <- colGrad[findInterval(subTable$density, colBins)]
    #pSize <- 3*(log10(subTable$size)/max(log10(subTable$size)))
    #points(xVals, yVals, pch = 16, cex = pSize, col = pCols)
    points(centromeres[i],cenY,pch=16,cex=3,bg="black")
    mtext(paste0("Chr", i),side=4,cex=1,adj=0.5,padj=-3)
    axis(1, at=seq(0,30,by=0.5),labels=xAchse,col.axis="black",las=0, cex.axis = 1.5, cex.lab = 1.5)
    if (i == 5) { mtext("position (Mb)",side=1,cex=1,adj=0.5,padj=4) }
    #if (i == "Chr2") { legend("right", title = "", legend = genotypes, col = unlist(genotypeCol[genotypes]), pch = 16, bty = 'n', cex =1.5, inset = 0.05, pt.cex = 1.5) }
  }
  dev.off()
  system(paste0("rsvg-convert -a -d 300 -p 300 ", file.path(rDir, outfileSVG)," > ", file.path(rDir, outfilePNG)))
  #system(paste("rm ", file.path(rDir, outfile), sep = ''))
}

###########################################################################################
####### Annotation related
###########################################################################################

#' Helper for \code{\link{f.get.genes2xyz.via.files}}
#'@param curFile current file to be read (key;entry|entryA|entryB)
#'@param singleCore use only one core
#'@return list[[key]] <- c(entryA, entryB)
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.get.genes2xyz.via.files.helper <- function(curFile, singleCore = FALSE) {
  if (!singleCore) {require("parallel")}
  temp <- scan(curFile, what = "character")
  temp <- unlist(strsplit(temp, split = ';', fixed =TRUE))
  out <- as.list(temp[seq(2, length(temp), by = 2)])
  names(out) <- temp[seq(1, length(temp), by = 2)]
  if (singleCore) {
    out <- lapply(out, function(x) unlist(strsplit(x, "|", fixed=TRUE)))
  } else {
    out <- mclapply(out, function(x) unlist(strsplit(x, "|", fixed=TRUE)), mc.cores =4)
  }
  return(out)
}

#' Get a genes2xyz list (from files prepared with makeMappingFiles.py)
#'@param gf2dmcFile AGI,feature;chrom_pos|chrom_pos (AT1G15050,exon;1_5182938|1_5182939)
#'@param g2gfFile AGI;AGI,feature|AGI,feature (AT1G06623;AT1G06623,five_flank|AT1G06623,miRNA_primary_transcript|AT1G06623,three_flank)
#'@param g2numFile AGI;number (AT1G06623;54)
#'@param singleCore use only one core
#'@return genes2xyz list
#'@seealso \code{\link{f.get.genes2xyz}} and \code{\link{f.get.transposons2xyz}}
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.get.genes2xyz.via.files <- function(gf2dmcFile, g2gfFile, g2numFile, singleCore = FALSE) {
  out <- list()
  out$dmc <- list()
  out$feature <- list()
  out$dmc <- f.get.genes2xyz.via.files.helper(gf2dmcFile)
  out$feature <- f.get.genes2xyz.via.files.helper(g2gfFile)
  out$numDMC <- as.numeric(unlist(f.get.genes2xyz.via.files.helper(g2numFile)))
  names(out$numDMC) <- names(f.get.genes2xyz.via.files.helper(g2numFile))
  return(out)
}

#' Get a genes2xyz list (from mappings)
#'@param dmcTable a data frame with the positions ("mapping" is required)
#'@param singleCore use only one core
#'@return genes2xyz list
#'@seealso \code{\link{f.get.genes2xyz.via.files}} and \code{\link{f.get.transposons2xyz}}
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.get.genes2xyz <- function(dmcTable, singleCore = FALSE) {
  if (!singleCore) {require("parallel")}
  out <- list()
  out$dmc <- list()
  out$feature <- list()
  ## extract genes that are DMC -> extract the features together with it and then create a list with AGI,feature->DMCs (use rownames)
  # extract gene,feature
  t.genes <- unique(grep("AT[[:alnum:]]{1}G[[:digit:]]{5}", unlist(strsplit(as.character(dmcTable$mapping), split = '|', fixed =TRUE)), value = TRUE))
  # extract the rownames of the DMCs
  if (singleCore) {
    genes2dmc <- lapply(t.genes, function(x) rownames(dmcTable)[grep(x, as.character(dmcTable$mapping))])
  } else {
    genes2dmc <- mclapply(t.genes, function(x) rownames(dmcTable)[grep(x, as.character(dmcTable$mapping))], mc.cores = 4)
  }
  # set the names in the list
  names(genes2dmc) <- t.genes
  ## get a list with gene->gene,feature
  #get the genes
  t.genes <- unique(grep("AT[[:alnum:]]{1}G[[:digit:]]{5}", unlist(strsplit(names(genes2dmc), split = ',', fixed =TRUE)), value = TRUE))
  # extract the gene,feature entries for each gene
  if (singleCore) {
    genes2gf <- lapply(t.genes, function(x) grep(x, names(genes2dmc), value = TRUE))
  } else {
    genes2gf <- mclapply(t.genes, function(x) grep(x, names(genes2dmc), value = TRUE), mc.cores = 4)
  }
  # set the names
  names(genes2gf) <- t.genes
  ## glue them together
  out$dmc <- genes2dmc
  out$feature <- genes2gf
  gfSum <- sapply(out$dmc, length)
  gSum <- sapply(out$feature, function(x) sum(gfSum[x]))
  out$numDMC <- gSum
  return(out)
}

#' Remove transposable element genes (TEGs) from genes2xyz lists
#'@param genes2xyz a list with gene-ID to sth mappings (see \code{\link{f.get.genes2xyz}}, \code{\link{f.get.genes2xyz.via.files}})
#'@param TEGtoTEmapFile tab separated file with two columns: TEG-AGI and TE-AGI
#'@return genes2xyz without TEGs (transposable element genes)
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.remove.tegs.from.genes2xyz <- function(genes2xyz, TEGtoTEmapFile) {
  teg2te <- read.table(TEGtoTEmapFile, header = FALSE, row.names = NULL, sep = '\t', stringsAsFactors = FALSE)
  tegs <- unique(teg2te[,1])
  out <- list(dmc = list(), feature = list())
  out$feature <- genes2xyz$feature[setdiff(names(genes2xyz$feature), tegs)]
  out$dmc <- genes2xyz$dmc[unlist(out$feature)]
  gfSum <- sapply(out$dmc, length)
  out$numDMC <- sapply(out$feature, function(x) sum(gfSum[x]))
  return(out)
}

#' Get a transposon2xyz list (from mappings)
#'@param dmcTable a data frame with the positions ("mapping" is required)
#'@param singleCore use only one core
#'@return transposon2xyz list
#'@seealso \code{\link{f.get.genes2xyz}} and \code{\link{f.get.genes2xyz.via.files}}
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.get.transposons2xyz <- function(dmcTable, singleCore = FALSE) {
  if (!singleCore) {require("parallel")}
  out <- list()
  out$dmc <- list()
  out$feature <- list()
  ## extract transposons that are DMC -> extract the features together with it and then create a list with AGI,feature->DMCs (use rownames)
  # extract gene,feature
  t.transposons <- unique(grep("AT[[:alnum:]]{1}TE[[:digit:]]{5}", unlist(strsplit(as.character(dmcTable$mapping), split = '|', fixed =TRUE)), value = TRUE))
  # extract the rownames of the DMCs
  if (singleCore) {
    transposons2dmc <- lapply(t.transposons, function(x) rownames(dmcTable)[grep(x, as.character(dmcTable$mapping))])
  } else {
    transposons2dmc <- mclapply(t.transposons, function(x) rownames(dmcTable)[grep(x, as.character(dmcTable$mapping))], mc.cores = 4)
  }
  # set the names in the list
  names(transposons2dmc) <- t.transposons
  ## get a list with gene->gene,feature
  #get the transposons
  t.transposons <- unique(grep("AT[[:alnum:]]{1}TE[[:digit:]]{5}", unlist(strsplit(names(transposons2dmc), split = ',', fixed =TRUE)), value = TRUE))
  # extract the gene,feature entries for each gene
  if (singleCore) {
    transposons2gf <- lapply(t.transposons, function(x) grep(x, names(transposons2dmc), value = TRUE))
  } else {
    transposons2gf <- mclapply(t.transposons, function(x) grep(x, names(transposons2dmc), value = TRUE), mc.cores = 4)
  }
  # set the names
  names(transposons2gf) <- t.transposons
  ## glue them together
  out$dmc <- transposons2dmc
  out$feature <- transposons2gf
  return(out)
}

###########################################################################################
####### Distance to annotation related
###########################################################################################

#' draw distance to an annotation feature (no random controls)
#'@param inDir directory where the distance files are located
#'@param rDir the directory in which the results will be stored
#'@param feature the annotation feature
#'@param smooth smoothen the line
#'@param outfilePrefix a prefix for the outfile name (default is "")
#'@param infilePrefix a prefix for the infile name (default is "")
#'@param summaryFunction function for summarization within bins
#'@return NULL
#'@note a distance file is named like this: <infilePrefix>_<feature>_<context>.txt
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.draw.distance.to.something <- function(inDir, rDir, feature, smooth = TRUE, outfilePrefix = "", infilePrefix = "", summaryFunction = sum) {
  svgOutfile <- file.path(rDir, paste0(outfilePrefix, "distance_to_", feature,".svg"))
  pngOutfile <- file.path(rDir, paste0(outfilePrefix, "distance_to_", feature,".png"))
  svg(svgOutfile, height = 10, width = 10)
  layout(matrix(1:4, nrow=2, byrow = TRUE))
  for (context in c("ALL", "CG", "CHG", "CHH")) {
    fileName <- paste0(infilePrefix, feature, "_", context, ".txt")
    distances <- read.table(file.path(inDir, fileName), sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
    if (smooth) {
      # bins on log scale - only works well with mean (considering that the bins are not of equal size on the linear scale) - well in a way it is not really true though
      distances$logDist <- log10(distances$dist+1)
      intervals <- seq(min(distances$logDist), max(distances$logDist), by = 0.1)
      bins <- split(distances, findInterval(distances$logDist, intervals, all.inside = FALSE))
      toUse <- as.numeric(names(bins))
      xVal <- intervals[toUse]
      yVal <- unlist(lapply(bins,function(x) ifelse(nrow(x) > 1, summaryFunction(x[,"occurence"]), x)))
      # bins on linear scale and then log (weird)
      # 				intervals <- seq(0, max(distances$dist), length.out = 100)
      # 				xVal <- log10(intervals+1)[1:99]
      # 				yVal <- split(distances, findInterval(distances$dist, intervals, all.inside = TRUE))
      # 				yVal <- unlist(lapply(yVal, function(x) ifelse(nrow(x) > 1, sum(x[,"occurence"]), x)))
    } else {
      xVal <- log10(distances$dist+1)
      yVal <- distances$occurence
    }
    f.histogram(xVal, yVal, xName = paste("distance to next", feature, "(log10())"), doNotShowSummary = TRUE, main = context, xlim=c(0,6))
  }
  dev.off()
  system(paste("rsvg-convert -a -d 300 -p 300 ", svgOutfile," > ", pngOutfile, sep = ''))
  return(NULL)
}

#' internal function for \code{\link{f.draw.distance.to.something.combined}}
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
f.draw.distance.to.something.get.random.mean.sd <- function(inDir, infilePrefix, feature, context, smooth = TRUE, summaryFunction = sum, numRep = 10, randomInDir = inDir) {
  # use unequal log steps (0.1 for 0 to 4 and then 0.2)
  #logSteps <- c(0, seq(1,2.8,by=0.2), seq(3,5,by=0.1))
  logSteps <- c(0, seq(1,5,by=0.5))
  temp <- list()
  maxX <- 0
  for (randomSetPrefix in paste("random", 1:numRep, sep = '_')) {
    cat(paste("processing ", randomSetPrefix, "\n", sep = ''))
    fileName <- paste0(randomSetPrefix, "_", infilePrefix, feature, "_", context, ".txt")
    distances <- read.table(file.path(randomInDir, fileName), sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
    #cat("limiting to 10'000\n")
    #distances <- subset(distances, dist < 1e4)
    if (smooth) {
      # bins on log scale - only works well with mean (considering that the bins are not of equal size on the linear scale) - well in a way it is not really true though
      distances$logDist <- log10(distances$dist+1)
      #intervals <- seq(min(distances$logDist), max(distances$logDist), by = 0.1)
      intervals <- logSteps
      bins <- split(distances, findInterval(distances$logDist, intervals, all.inside = TRUE))
      toUse <- as.numeric(names(bins))
      xVal <- intervals[toUse]
      yVal <- unlist(lapply(bins,function(x) ifelse(nrow(x) > 1, summaryFunction(x[,"occurence"]), x)))
      # bins on linear scale and then log (weird)
      # 				intervals <- seq(0, max(distances$dist), length.out = 100)
      # 				xVal <- log10(intervals+1)[1:99]
      # 				yVal <- split(distances, findInterval(distances$dist, intervals, all.inside = TRUE))
      # 				yVal <- unlist(lapply(yVal, function(x) ifelse(nrow(x) > 1, sum(x[,"occurence"]), x)))
    } else {
      xVal <- log10(distances$dist+1)
      yVal <- distances$occurence
    }
    maxX <- ifelse(max(xVal) > maxX, max(xVal), maxX)
    # scale all to relative occurences because we just have the random data for the total
    temp[[randomSetPrefix]] <- cbind(xVal, yVal/sum(yVal))
  }
  #rn <- as.character(seq(0, maxX, by = 0.1))
  rn <- as.character(logSteps)
  mat <- matrix(0, nrow = length(rn), ncol = numRep); colnames(mat) <- names(temp); rownames(mat) <- rn
  for (randomSetPrefix in names(temp)) {
    mat[as.character(temp[[randomSetPrefix]][,1]),randomSetPrefix] <- temp[[randomSetPrefix]][,2]
  }
  out <- cbind(apply(mat, 1, mean), apply(mat, 1, sd))
  return(out)
}

#' draw distance to an annotation feature (with random controls); line graphs
#'@param inDir directory where the distance files are located
#'@param rDir the directory in which the results will be stored
#'@param feature the annotation feature
#'@param smooth smoothen the line
#'@param outfilePrefix a prefix for the outfile name (default is "")
#'@param infilePrefix a prefix for the infile name (default is "")
#'@param summaryFunction function for summarization within bins
#'@param numRep number of random samples
#'@param randomInDir directory where the RANDOM distance files are located
#'@return NULL
#'@note a distance file is named like this: <infilePrefix>_<feature>_<context>[_<gain/loss>].txt
#'the randomized files must be: random_<number>_<infilePrefix>_<feature>_<context>[_<gain/loss>].txt
#'[_<gain/loss>] is optional (without it, it will assume "totalC")
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.draw.distance.to.something.combined <- function(inDir, rDir, feature, smooth = TRUE, outfilePrefix = "", infilePrefix = "",
                                                  summaryFunction = sum, numRep = 10, randomInDir = inDir) {
  # use unequal log steps (0.1 for 0 to 4 and then 0.2)
  #logSteps <- c(0, seq(1,2.8,by=0.2), seq(3,5,by=0.1))
  logSteps <- c(0, seq(1,5,by=0.5))
  contextColors <- data.frame(
    ALL = c("black", "green", "red"),
    CG = brewer.pal(3,"YlOrBr"),
    CHG = brewer.pal(3,"PuBuGn"),
    CHH = brewer.pal(3,"BuGn"),
    row.names = c("totalC", "gain", "loss"), stringsAsFactors = FALSE
  )
  cat("note that this throws an error if you use smooth = FALSE\n")
  svgOutfile <- file.path(rDir, paste0(outfilePrefix, "distance_to_", feature,".svg"))
  pngOutfile <- file.path(rDir, paste0(outfilePrefix, "distance_to_", feature,".png"))
  svg(svgOutfile, height = 10, width = 10)
  layout(matrix(1:4, nrow=2, byrow = TRUE))
  for (context in c("ALL", "CG", "CHG", "CHH")) {
    xVals <- list()
    yVals <- list()
    for (gainLoss in c("totalC", "gain", "loss")) {
      if (gainLoss == "totalC") {
        fileName <- file.path(inDir, paste0(infilePrefix, feature, "_", context, ".txt"))
      } else {
        fileName <- file.path(inDir, paste0(infilePrefix, feature, "_", context, "_", gainLoss, ".txt"))
      }
      if (!file.exists(fileName)) { cat("skipping:", fileName, "\n"); next }
      distances <- read.table(fileName, sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
      #cat("limiting to 10'000\n")
      #distances <- subset(distances, dist < 1e4)
      if (smooth) {
        # bins on log scale - only works well with mean (considering that the bins are not of equal size on the linear scale) - well in a way it is not really true though
        distances$logDist <- log10(distances$dist+1)
        #intervals <- seq(min(distances$logDist), max(distances$logDist), by = 0.1)
        intervals <- logSteps
        bins <- split(distances, findInterval(distances$logDist, intervals, all.inside = TRUE))
        toUse <- as.numeric(names(bins))
        xVals[[gainLoss]] <- intervals[toUse]
        yVals[[gainLoss]] <- unlist(lapply(bins,function(x) ifelse(nrow(x) > 1, summaryFunction(x[,"occurence"]), x)))
        # bins on linear scale and then log (weird)
        # 				intervals <- seq(0, max(distances$dist), length.out = 100)
        # 				xVal <- log10(intervals+1)[1:99]
        # 				yVal <- split(distances, findInterval(distances$dist, intervals, all.inside = TRUE))
        # 				yVal <- unlist(lapply(yVal, function(x) ifelse(nrow(x) > 1, sum(x[,"occurence"]), x)))
        # use density function
      } else {
        xVals[[gainLoss]] <- log10(distances$dist+1)
        yVals[[gainLoss]] <- distances$occurence
      }
      # scale all to relative occurences because we just have the random data for the total
      yVals[[gainLoss]] <- yVals[[gainLoss]]/sum(yVals[[gainLoss]])
    }
    randomData <- f.draw.distance.to.something.get.random.mean.sd(inDir, infilePrefix, feature, context, smooth, summaryFunction, numRep, randomInDir)
    obsMax <- max(unlist(yVals))
    randMax <- max(randomData[,1]+1.96*abs(randomData[,2]))
    maxYval <- ifelse(obsMax > randMax, obsMax, randMax)
    maxYval <- ifelse(maxYval > 0.45, maxYval, 0.45)
    f.histogram(xVals$totalC, yVals$totalC, xName = paste("distance to next", feature, "(log10())"), doNotShowSummary = TRUE,
                main = context, xlim=c(0,5), ylim = c(0, maxYval), lwd = 2, col = contextColors["totalC",context])
    xValRand <- as.numeric(rownames(randomData))
    lines(xValRand, randomData[,1]+1.96*randomData[,2], lty="dashed", col="black", lwd=1)
    lines(xValRand, randomData[,1], lty="solid", col="black", lwd=2)
    lines(xValRand, randomData[,1]-1.96*randomData[,2], lty="dashed", col="black", lwd=1)
    xPolygon <- c(xValRand, rev(xValRand))
    yPolygon <- c(randomData[,1]+1.96*randomData[,2],rev(randomData[,1]-1.96*randomData[,2]))
    polygon(xPolygon, yPolygon, col="#00009925", border=NA)
    if ("gain" %in% names(xVals)) { lines(xVals$gain, yVals$gain, lty="solid", col = contextColors["gain",context], lwd=2) }
    if ("loss" %in% names(xVals)) { lines(xVals$loss, yVals$loss, lty="solid", col = contextColors["loss",context], lwd=2) }
  }
  dev.off()
  system(paste("rsvg-convert -a -d 300 -p 300 ", svgOutfile," > ", pngOutfile, sep = ''))
}

#' draw distance to an annotation feature (with random control); boxplots
#'@param inDir directory where the distance files are located
#'@param rDir the directory in which the results will be stored
#'@param feature the annotation feature
#'@param outfilePrefix a prefix for the outfile name (default is "boxplots_")
#'@param infilePrefix a prefix for the infile name (default is "")
#'@param randomInDir directory where the RANDOM distance files are located
#'@return NULL
#'@note a distance file is named like this: <infilePrefix>_<feature>_<context>[_<gain/loss>].txt
#'the randomized files must be: random_<number>_<infilePrefix>_<feature>_<context>[_<gain/loss>].txt
#'[_<gain/loss>] is optional (without it, it will assume "totalC")
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.draw.distance.to.something.boxplot <- function(inDir, rDir, feature, outfilePrefix = "boxplots_", infilePrefix = "", randomInDir = inDir) {
  contextColors <- data.frame(
    ALL = c("black", "green", "red"),
    CG = brewer.pal(3,"YlOrBr"),
    CHG = brewer.pal(3,"PuBuGn"),
    CHH = brewer.pal(3,"BuGn"),
    row.names = c("totalC", "gain", "loss"), stringsAsFactors = FALSE
  )
  boxplotColors <- c()
  svgOutfile <- file.path(rDir, paste(outfilePrefix, "distance_to_", feature,".svg", sep = ''))
  pngOutfile <- file.path(rDir, paste(outfilePrefix, "distance_to_", feature,".png", sep = ''))
  svg(svgOutfile, height = 5, width = 5)
  toPlot <- list()
  for (context in c("CG", "CHG", "CHH")) {
    fileName <- paste("random_1_", infilePrefix, feature, "_", context, ".txt", sep = '')
    distances <- read.table(file.path(randomInDir, fileName), sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
    revDat <- log10(f.reverse.histogram.data.matrix(as.matrix(distances))+1)
    groupName <- paste("random", context, sep = '_')
    toPlot[[groupName]] <- data.frame(distance = revDat, group = rep(groupName, length(revDat)))
    boxplotColors <- c(boxplotColors, "gray30")
    for (gainLoss in c("totalC", "gain", "loss")) {
      if (gainLoss == "totalC") {
        groupName <- paste(feature, context, sep = '_')
        fileName <- file.path(inDir, paste0(infilePrefix, feature, "_", context, ".txt"))
      } else {
        groupName <- paste(feature, context, gainLoss, sep = '_')
        fileName <- file.path(inDir, paste0(infilePrefix, feature, "_", context, "_", gainLoss, ".txt"))
      }
      if (!file.exists(fileName)) { cat("skipping:", fileName, "\n"); next}
      distances <- read.table(fileName, sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
      revDat <- log10(f.reverse.histogram.data.matrix(as.matrix(distances))+1)
      toPlot[[groupName]] <- data.frame(distance = revDat, group = rep(groupName, length(revDat)))
      boxplotColors <- c(boxplotColors, contextColors[gainLoss,context])
    }
  }
  toPlot <- do.call("rbind", toPlot)
  boxplot(distance ~ group, data = toPlot, col = boxplotColors, pch = 16, cex = 0.5, las = 1, xaxs = "r", yaxs = "r", bty = "n", ylim = c(0, 5))
  dev.off()
  system(paste("rsvg-convert -a -d 300 -p 300 ", svgOutfile," > ", pngOutfile, sep = ''))
  return(NULL)
}

#' test distance to an annotation feature - random data - t-test
#'@param inDir directory where the distance files are located
#'@param feature the annotation feature
#'@param infilePrefix a prefix for the infile name (default is "")
#'@param numRep number of random samples
#'@param randomInDir directory where the RANDOM distance files are located
#'@return a table with the results
#'@note a distance file is named like this: <infilePrefix>_<feature>_<context>[_<gain/loss>].txt
#'the randomized files must be: random_<number>_<infilePrefix>_<feature>_<context>[_<gain/loss>].txt
#'[_<gain/loss>] is optional (without it, it will assume "totalC")
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.test.distance.to.something <- function(inDir, feature, infilePrefix = "", numRep = 10, randomInDir = inDir) {
  out <- list()
  for (context in c("ALL", "CG", "CHG", "CHH")) {
    for (gainLoss in c("totalC", "gain", "loss")) {
      if (gainLoss == "totalC") {
        groupName <- paste(feature, context, sep = '_')
        fileName <- file.path(inDir, paste0(infilePrefix, feature, "_", context, ".txt"))
      } else {
        groupName <- paste(feature, context, gainLoss, sep = '_')
        fileName <- file.path(inDir, paste0(infilePrefix, feature, "_", context, "_", gainLoss, ".txt"))
      }
      if (!file.exists(fileName)) { cat("skipping:", fileName, "\n"); next }
      observedDistanceData <- read.table(fileName, sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
      observedDistances <- log10(f.reverse.histogram.data.matrix(as.matrix(observedDistanceData))+1)
      for (randomSetPrefix in paste("random", 1:numRep, sep = '_')) {
        fileName <- paste0(randomSetPrefix, "_", infilePrefix, feature, "_", context, ".txt")
        randomDistanceData <- read.table(file.path(randomInDir, fileName), sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
        randomDistances <- log10(f.reverse.histogram.data.matrix(as.matrix(randomDistanceData))+1)
        temp <- t.test(observedDistances, randomDistances, alternative = "two.sided")
        out[[paste(context, gainLoss, randomSetPrefix, sep = '_')]] <- c(temp$estimate,temp$p.value)
      }
    }
  }
  rn <- names(out)
  out <- do.call("rbind", out)
  rownames(out) <- rn
  colnames(out) <- c("observedMean", "randomMean", "pValue")
  return(out)
}

#' test distance to an annotation feature - random data - empirical (needs high numRep)
#'@param inDir directory where the distance files are located
#'@param feature the annotation feature
#'@param infilePrefix a prefix for the infile name (default is "")
#'@param numRep number of random samples
#'@param randomInDir directory where the RANDOM distance files are located
#'@return a table with the results
#'@note a distance file is named like this: <infilePrefix>_<feature>_<context>[_<gain/loss>].txt
#'the randomized files must be: random_<number>_<infilePrefix>_<feature>_<context>[_<gain/loss>].txt
#'[_<gain/loss>] is optional (without it, it will assume "totalC")
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.test.distance.to.something.empirical <- function(inDir, feature, infilePrefix = "", numRep = 500, randomInDir = inDir) {
  cat("context == ALL will be skipped\n")
  out <- list()
  for (context in c("CG", "CHG", "CHH")) {
    for (gainLoss in c("totalC", "gain", "loss")) {
      if (gainLoss == "totalC") {
        groupName <- paste(feature, context, sep = '_')
        fileName <- file.path(inDir, paste0(infilePrefix, feature, "_", context, ".txt"))
      } else {
        groupName <- paste(feature, context, gainLoss, sep = '_')
        fileName <- file.path(inDir, paste0(infilePrefix, feature, "_", context, "_", gainLoss, ".txt"))
      }
      if (!file.exists(fileName)) { cat("skipping:", fileName, "\n"); next }
      observedDistanceData <- read.table(file.path(inDir, fileName), sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
      observedDistances <- log10(f.reverse.histogram.data.matrix(as.matrix(observedDistanceData))+1)
      numC <- length(observedDistances)
      aveObs <- mean(observedDistances)
      aveRan <- c()
      for (randomSetPrefix in paste("random", 1:numRep, sep = '_')) {
        fileName <- paste0(randomSetPrefix, "_", infilePrefix, feature, "_", context, ".txt")
        randomDistanceData <- read.table(file.path(randomInDir, fileName), sep = "\t", quote = "", header = FALSE, row.names = NULL, col.names = c("dist", "occurence"))
        randomDistances <- log10(sample(f.reverse.histogram.data.matrix(as.matrix(randomDistanceData)), numC)+1)
        aveRan <- c(aveRan, mean(randomDistances))
      }
      out[[paste(genotype, "_", context, "_", gainLoss, sep = '')]] <- c(aveObs, mean(aveRan), sd(aveRan), sum(aveRan < aveObs)/numRep)
      cat(paste(sum(aveRan < aveObs)/numRep, "\n", sep = ''))
    }
  }
  rn <- names(out)
  out <- do.call("rbind", out)
  rownames(out) <- rn
  colnames(out) <- c("observedMean", "randomMean", "randomSd", "pValue")
  return(out)
}

###########################################################################################
####### Number of DMCs - a plot
###########################################################################################

#' Draw number of DMCs per gene
#'@param genes2xyz a list with gene-ID to sth mappings (see \code{\link{f.get.genes2xyz}}, \code{\link{f.get.genes2xyz.via.files}} and \code{\link{f.get.transposons2xyz}})
#'@param rDir the directory in which the results will be stored
#'@param filePrefix a prefix for the file name (default is "")
#'@param minimalNumber
#'@return NULL
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.plot.number.DMCs.per.gene <- function(genes2xyz, rDir, filePrefix = "", minimalNumber = 1) {
  if ("numDMC" %in% names(genes2xyz)) {
    gSum <- genes2xyz$numDMC
  } else {
    # for each gf, get the number of dmcs
    gfSum <- sapply(genes2xyz$dmc, length)
    # sum up for each gene
    gSum <- sapply(genes2xyz$feature, function(x) sum(gfSum[x]))
  }
  svgOutfile <- file.path(rDir, paste(filePrefix, "DMCs_per_gene.svg", sep = ''))
  pngOutfile <- file.path(rDir, paste(filePrefix, "DMCs_per_gene.png", sep = ''))
  svg(svgOutfile, height = 5, width = 5)
  f.histogram(log2(gSum[gSum >= minimalNumber]), xName = "number of DMCs per gene (log2(x+1))", doNotShowSummary = TRUE, main = "", xlim = c(log2(minimalNumber), 12))
  dev.off()
  system(paste("rsvg-convert -a -d 300 -p 300 ", svgOutfile," > ", pngOutfile, sep = ''))
  write.table(cbind(names(gSum), gSum), file.path(rDir, paste0(filePrefix, "_DMCs_per_gene.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  return(NULL)
}

###########################################################################################
####### Distances between Cs related
###########################################################################################

#' Calculate distances between Cs
#'@param dmcTable a data frame with the positions ("chrom" and "pos" are required)
#'@return a vector with all distances
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.get.distances.between <- function(dmcTable) {
  out <- c()
  for (chr in unique(dmcTable$chrom)) {
    sub <- subset(dmcTable, chrom == chr)
    num <- nrow(sub)
    if (num >= 2) { out <- c(out, log10(sub$pos[2:num] - sub$pos[1:(num-1)])) }
  }
  return(out)
}

#' Draw distances between Cs
#'@param dmcTable a data frame with the positions ("chrom" and "pos" are required)
#'@param rDir the directory in which the results will be stored
#'@param filePrefix a prefix for the file name (default is "")
#'@return NULL
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.draw.distance.between.DMCs <- function(dmcTable, rDir, filePrefix = "") {
  myDist <- f.get.distances.between(dmcTable)
  svgOutfile <- file.path(rDir, paste(filePrefix, "distances_between_DMCs.svg", sep = ''))
  pngOutfile <- file.path(rDir, paste(filePrefix, "distances_between_DMCs.png", sep = ''))
  svg(svgOutfile, height = 5, width = 5)
  f.histogram(myDist, xName = "distance to next DMC (log10())", doNotShowSummary = TRUE, main = "", xlim=c(0,6))
  dev.off()
  system(paste("rsvg-convert -a -d 300 -p 300 ", svgOutfile," > ", pngOutfile, sep = ''))
  return(NULL)
}

#' Draw distances between RANDOM Cs
#'@param dmcTable a data frame with the positions ("chrom" and "pos" are required)
#'@param allPositions a data frame with all possible position (columns "chrom" and "pos" are required)
#'@param rDir the directory in which the results will be stored
#'@param numRep number of random samples
#'@return NULL
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.draw.distance.between.randomCs <- function(dmcTable, allPositions, rDir, numRep) {
  randomDMC <- dmcTable[,c("chrom", "pos")]
  for (randomSetPrefix in paste("random", 1:numRep, sep = '_')) { ## NOTE THAT THIS HERE IS DEFINITELY NOT OPTIMIZED FOR LARGE NUMBER OF SETS
    cat(paste("processing ", randomSetPrefix, "\n", sep = ''))
    for (chr in unique(dmcTable$chrom)) {
      randomDMC$pos[randomDMC$chrom == chr] <- sort(sample(allPositions$pos[allPositions$chrom == chr], sum(dmcTable$chrom == chr)))
    }
    cat(paste("sampled ", randomSetPrefix, "\n", sep = ''))
    f.draw.distance.between.DMCs(randomDMC, rDir, randomSetPrefix)
  }
}

#' get random mean and sd for drawing random and observed distances
#'@param dmcTable a data frame with the positions ("chrom" and "pos" are required)
#'@param allPositions a data frame with all possible position (columns "chrom" and "pos" are required)
#'@param numRep number of random samples
#'@return matrix with random means and standard deviations
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.draw.distance.between.DMCs.combined.get.random.mean.sd <- function(dmcTable, allPositions, numRep) {
  chromCounts <- table(allPositions$chrom)
  chromsToKeep <- names(chromCounts)[chromCounts > 1]
  dmcTable <- subset(dmcTable, chrom %in% chromsToKeep)
  allPositions <- subset(allPositions, chrom %in% chromsToKeep)
  randomDMC <- dmcTable[,c("chrom", "pos")]
  randomDensities <- list()
  for (randomSetPrefix in paste("random", 1:numRep, sep = '_')) {
    cat(paste("processing ", randomSetPrefix, "\n", sep = ''))
    for (chr in unique(randomDMC$chrom)) {
      randomDMC$pos[randomDMC$chrom == chr] <- sort(sample(allPositions$pos[allPositions$chrom == chr], sum(randomDMC$chrom == chr)))
    }
    DMCdistances <- f.get.distances.between(randomDMC)
    randomDensities[[randomSetPrefix]] <- density(DMCdistances, from = 0, to = 6, n = 50)$y
    cat(paste("sampled ", randomSetPrefix, "\n", sep = ''))
  }
  randomDensities <- do.call("cbind", randomDensities)
  out <- cbind(apply(randomDensities, 1, mean), apply(randomDensities, 1, sd))
  rownames(out) <- as.character(seq(0, 6, length.out = nrow(out)))
  return(out)
}

#' draw observed and random distances between positions
#'@param dmcTable a data frame with the positions ("chrom" and "pos" are required)
#'@param allPositions a data frame with all possible position (columns "chrom" and "pos" are required)
#'@param rDir the directory in which the results will be stored
#'@param numRep number of random samples
#'@param filePrefix a prefix for the file name (default is "")
#'@return NULL
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.draw.distance.between.DMCs.combined <- function(dmcTable, allPositions, rDir, numRep, filePrefix = "") {
  DMCdistances <- f.get.distances.between(dmcTable)
  randomData <- f.draw.distance.between.DMCs.combined.get.random.mean.sd(dmcTable, allPositions, numRep)
  svgOutfile <- file.path(rDir, paste(filePrefix, "distances_between_DMCs.svg", sep = ''))
  pngOutfile <- file.path(rDir, paste(filePrefix, "distances_between_DMCs.png", sep = ''))
  svg(svgOutfile, height = 5, width = 5)
  plot(NULL, type = "n", xlim = c(0, 6), ylim = c(0, 1), xlab = "distance to next DMC (log10())", ylab = "density", main = "", bty = "n", las = 1)
  obsDensities <- density(DMCdistances, from = 0, to = 6, n = 50)
  lines(obsDensities$x, obsDensities$y, type = 'l', lwd = 2)
  xValRand <- as.numeric(rownames(randomData))
  lines(xValRand, randomData[,1]+1.96*randomData[,2], lty="dashed", col="black", lwd=1)
  lines(xValRand, randomData[,1], lty="solid", col="black", lwd=2)
  lines(xValRand, randomData[,1]-1.96*randomData[,2], lty="dashed", col="black", lwd=1)
  xPolygon <- c(xValRand, rev(xValRand))
  yPolygon <- c(randomData[,1]+1.96*randomData[,2],rev(randomData[,1]-1.96*randomData[,2]))
  polygon(xPolygon, yPolygon, col="#00009925", border=NA)
  dev.off()
  system(paste("rsvg-convert -a -d 300 -p 300 ", svgOutfile," > ", pngOutfile, sep = ''))
  return(NULL)
}

###########################################################################################
####### Region related
###########################################################################################

#' Create regions (e.g. DMRs) from candidate positions.
#'@param dmcTable a data frame with the positions (required are the columns "chrom" and "pos")
#'@param minNum minimal number of positions within a region ("chrom", "start", "end", "size", "density", "numbers")
#'@param maxDist maximal distance between two neighboring positions to be joined into a region
#'@return a data frame with the regions
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.create.regions <- function(dmcTable, minNum, maxDist) {
  # could be optimized if necessary
  starts <- c()
  ends <- c()
  numbers <- c()
  chromosomes <- unique(dmcTable$chrom)
  regionCounts <- rep(0, length(chromosomes))
  names(regionCounts) <- chromosomes
  for (chr in chromosomes) {
    sub <- subset(dmcTable, chrom == chr)
    numPos <- nrow(sub)
    if (numPos < minNum) { next }
    curStart <- sub$pos[1]
    curEnd <- sub$pos[1]
    curNum <- 1
    for (i in 2:numPos) {
      nextPos <- sub$pos[i]
      if ((nextPos-curEnd) <= maxDist) {
        curEnd <- nextPos
        curNum <- curNum + 1
      } else {
        if (curNum >= minNum) { 
          starts <- c(starts, curStart)
          ends <- c(ends, curEnd)
          numbers <- c(numbers, curNum)
          regionCounts[chr] <- regionCounts[chr]+1
        }
        curStart <- nextPos
        curEnd <- nextPos
        curNum <- 1
      }
    }
    if (curNum >= minNum) { 
      starts <- c(starts, curStart)
      ends <- c(ends, curEnd)
      numbers <- c(numbers, curNum)
      regionCounts[chr] <- regionCounts[chr]+1
    }
  }
  size <- ends-starts
  out <- data.frame(chrom = rep(chromosomes, times = regionCounts), start = starts, end = ends, size = size, density = numbers/size, count = numbers, stringsAsFactors = FALSE)
  return(out)
}

#' Extract positions which are located in a defined region.
#'@param dmcTable a data frame with the positions (required are the columns "chrom" and "pos")
#'@param regionTable a data frame with the regions (required are the columns "chrom", "start", "end")
#'@return a boolean vector of the length nrow(dmcTable) [TRUE means the position is within]
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.extract.Cs.within.regions <- function(dmcTable, regionTable) {
  temp <- list()
  for (chr in unique(regionTable$chrom)) {
    subDMC <- subset(dmcTable, chrom == chr)
    subREG <- subset(regionTable, chrom == chr)
    temp[[chr]] <- unlist(apply(subREG[,c("start","end")], 1, function(x) rownames(subDMC)[subDMC$pos > x[1] & subDMC$pos < x[2]]))
  }
  out <- unlist(temp)
  return(out)
}

#' Join overlapping regions
#'@param data a data frame ("chrom", "start", "end" are requires)
#'@param otherCols a vector with column names to be summarized as well (using summaryFunction)
#'@param summaryFunction function to summarize the otherCols
#'@param gap allow a gap between regions
#'@return merged regions
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.join.overlapping.regions <- function(data, otherCols = c(), summaryFunction = mean, gap = 0) {
  data <- data[with(data, order(chrom, start, end)),]
  data$regNum <- 1
  offset <- 1
  for (chr in unique(data$chrom)) {
    subData <- subset(data, chrom == chr)
    subData$regNum[1] <- offset
    if (nrow(subData) > 1) {
      temp <- (subData$start[2:nrow(subData)] - subData$end[1:(nrow(subData)-1)]) > gap
      subData$regNum[2:nrow(subData)] <- offset + cumsum(temp)
    } 
    data[data$chrom == chr, ] <- subData
    offset <- max(data$regNum) + 1
  }
  out <- data.frame(
    chrom = aggregate(data$chrom, by=list(region=data$regNum), unique, simplify = TRUE)$x,
    start = aggregate(data$start, by=list(region=data$regNum), min, simplify = TRUE)$x,
    end = aggregate(data$end, by=list(region=data$regNum), max, simplify = TRUE)$x,
    count = aggregate(data$chrom, by=list(region=data$regNum), length, simplify = TRUE)$x,
    stringsAsFactors = FALSE
  )
  for (cn in otherCols) {
    out[[cn]] = aggregate(data[[cn]], by=list(region=data$regNum), summaryFunction, simplify = TRUE)$x
  }
  return(out)
}

#' Create a series of regions with different parameters and draw a heatmap as overview.
#'@param dmcTable a data frame with the positions (columns "chrom" and "pos" are required)
#'@param rDir the directory in which the results will be stored
#'@param filePrefix a prefix for the file name (default is "")
#'@param bySpacerNucleotides the minimal number of DMCs per region (a vector with several values to test)
#'@param bySpacerDistances the maximal distance between two neighboring DMCs (a vector with several values to test)
#'@param useLog draw log2(x+1) instead of counts
#'@return a list with some matrices (number of regions, number DMCs within regions, DMC density)
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.analyze.regions <- function(dmcTable, rDir, filePrefix = "", bySpacerNucleotides = seq(5, 100, by = 5), bySpacerDistances = seq(10, 200, by = 10), useLog = FALSE) {
  dataTypes <- c("numberRegions", "withinRegions", "density")
  xChars <- as.character(bySpacerNucleotides)
  yChars <- as.character(bySpacerDistances)
  out <- list()
  for (dt in dataTypes) {
    out[[dt]] <- matrix(0, nrow = length(bySpacerNucleotides), ncol = length(bySpacerDistances), dimnames = list(xChars, yChars))
  }
  for (minNuc in bySpacerNucleotides) {
    cat("processing", minNuc, "...\n")
    for (maxDist in bySpacerDistances) {
      temp <- f.create.regions(dmcTable, minNuc, maxDist)
      out$numberRegions[as.character(minNuc), as.character(maxDist)] <- nrow(temp)
      out$withinRegions[as.character(minNuc), as.character(maxDist)] <- sum(temp$count)
      out$density[as.character(minNuc), as.character(maxDist)] <- round(sum(temp$count)/nrow(temp), 0)
    }
  }
  svgOutfile <- file.path(rDir, paste(filePrefix, "region_parameters_and_results.svg", sep = ''))
  pngOutfile <- file.path(rDir, paste(filePrefix, "region_parameters_and_results.png", sep = ''))
  svg(svgOutfile, height = 15, width = 10)
  layout(matrix(1:3, ncol = 1))
  xPos <- 1:length(bySpacerNucleotides)
  yPos <- 1:length(bySpacerDistances)
  xLabs <- "minimal number of DMCs"
  yLabs <- "maximal distance between two DMCs"
  mainLabs <- c("number of regions", "number of DMCs within regions", "number of DMCs per region"); names(mainLabs) <- dataTypes
  for (dt in dataTypes) {
    toPlot <- out[[dt]]
    if (useLog & dt %in% c("numberRegions", "withinRegions")) {
      modPlot <- log2(toPlot + 1)
    } else {
      modPlot <- toPlot
    }
    image(xPos, yPos, modPlot, xlab = xLabs, ylab = yLabs, main = mainLabs[dt], yaxt = "n", xaxt = "n")
    axis(1, at = xPos, labels = xChars, outer = FALSE)
    axis(2, at = yPos, labels = yChars, outer = FALSE)
    text(rep(xPos, each = length(yPos)), yPos, t(toPlot), cex = 1.7)
  }
  dev.off()
  system(paste("rsvg-convert -a -d 300 -p 300 ", svgOutfile," > ", pngOutfile, sep = ''))
  return(out)
}

#' Create a series of RANDOM regions with different parameters and draw a heatmap as overview.
#'@param dmcTable a data frame with the significant positions (columns "chrom" and "pos" are required)
#'@param allPositions a data frame with all possible position (columns "chrom" and "pos" are required)
#'@param rDir the directory in which the results will be stored
#'@param numRep number of random samples
#'@return NULL
#'@note 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch}.
#'@export
f.analyze.random.regions <- function(dmcTable, allPositions, rDir, numRep) {
  randomDMC <- dmcTable[,c("chrom", "pos")]
  for (randomSetPrefix in paste("random", 1:numRep, sep = '_')) { ## NOTE THAT THIS HERE IS DEFINITELY NOT OPTIMIZED FOR LARGE NUMBER OF SETS
    cat(paste("processing ", randomSetPrefix, "\n", sep = ''))
    for (chr in unique(dmcTable$chrom)) {
      randomDMC$pos[randomDMC$chrom == chr] <- sample(allPositions$pos[allPositions$chrom == chr], sum(randomDMC$chrom == chr))
    }
    cat(paste("sampled ", randomSetPrefix, "\n", sep = ''))
    temp <- f.analyze.regions(randomDMC, rDir, randomSetPrefix)
  }
  return(NULL)
}

