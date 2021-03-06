#!/usr/bin/env Rscript
library("stats")
library("parallel")

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
myPath <- as.character(myarg[argPos+1])
numberOfCores <- as.numeric(myarg[argPos+2])
useDfCor <- as.character(myarg[argPos+3])

## print a message function
f.print.message <- function(x) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste0("=== ", x,"\n")) }

## specify the model used for the analysis of an individual cytosine
# the example assumes that you have two factors with two levels each (same example as in MethAn/MethAnPre/mergeBedGraphs.py):
# sampleID,group,organ,sex
# lm1,liver_male,liver,male
# lm2,liver_male,liver,male
# lm3,liver_male,liver,male
# lf1,liver_female,liver,female
# lf2,liver_female,liver,female
# lf3,liver_female,liver,female
# hm1,heart_male,heart,male
# hm2,heart_male,heart,male
# hm3,heart_male,heart,male
# hf1,heart_female,heart,female
# hf2,heart_female,heart,female
# hf3,heart_female,heart,female
# using mergeBedGraphs.py and the example given there, the columns of your input for this script would be:
# chrom, pos, ctxt, total coverage, percent methylation, sample, group, organ, sex
# the full model would be: percentMethylation ~ organ + sex + organ:sex
# below there are two functions implementing this model
# the first is the basic implementation with lm
# the second implementation is experimental and it tries to correct the df if fitted values are 0 or 100, it's entirely inspired by:
# Lun, ATL, and Smyth, GK (2017). No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data. Statistical Applications in Genetics and Molecular Biology 16.
## IMPORTANT:
# there is no check/rescue for incomplete data (i.e., the function assumes that there are >1 non-NAs for all groups, which is true if you use mergeBedGraphs.py)
# the function with the df-correction assumes that all explanatory variables are factors and NOT CONTINUOUS VARIABLES
# (df-correction can be done with continuous factors as well, but it does not work with the given example)
basicModel <- function(x){
  mod <- as.matrix(anova(lm(terms(pcentmeth~organ+sex+organ:sex, keep.order=TRUE), data = x)))
  out <- c(x[1,1:3], nrow(x), mean(x[,4], na.rm = TRUE), sd(x[,4], na.rm = TRUE), mean(x[,5], na.rm = TRUE), sd(x[,5], na.rm = TRUE),
           mod[1:3, "F value"], mod[1:3, "Pr(>F)"], mod["Residuals","Mean Sq"])
  return(out)
}
dfCorrectedModel <- function(x){
  fit <- lm(terms(pcentmeth~organ+sex+organ:sex, keep.order = TRUE),data = x)
  mod <- as.matrix(anova(fit))
  fitCloseToZero <- (fit$fitted.values < 0.0001)
  fitCloseToHundret <- (fit$fitted.values > 99.9999)
  fitCloseToZeroOrHundret <- fitCloseToZero | fitCloseToHundret
  addOne <- as.numeric((sum(fitCloseToZero)>0)&(sum(fitCloseToHundret)>0)) # add one df again if some are close to 0 and others are close to 100
  standardDF <- mod["Residuals", "Df"]
  temp <- fit$model[(!fitCloseToZeroOrHundret),]
  numCoeff <- length(unique(apply(temp[,2:ncol(temp)], 1, function(x) paste(x, collapse = '_'))))
  adjustedDF <- nrow(x) - sum(fitCloseToZeroOrHundret) + addOne - numCoeff
  if (standardDF != adjustedDF) {
    adjustedMS <- mod["Residuals", "Sum Sq"]/adjustedDF
    mod[1:3, "F value"] <- mod[1:3, "Mean Sq"]/adjustedMS
    mod[1:3, "Pr(>F)"] <- 1-pf(mod[1:3, "F value"], mod[1:3, "Df"], adjustedDF)
  } else {
    adjustedMS <- mod["Residuals", "Mean Sq"]
  }
  out <- c(x[1,1:3], nrow(x), mean(x[,4], na.rm = TRUE), sd(x[,4], na.rm = TRUE), mean(x[,5], na.rm = TRUE), sd(x[,5], na.rm = TRUE),
           mod[1:3, "F value"], mod[1:3, "Pr(>F)"], adjustedMS)
  return(out)
}

# this are the entries returned by the functions above
modelOutputNames <- c("chrom", "pos", "context", "numSamples", "meanCov", "sdCov", "meanPcent", "sdPcent", # fixed
                      paste0(rep(c("F_", "P_"), each = 3), c("organ", "sex", "organXsex")), "MS_res") # except for MS_res, this depends on the model

# specify the headers in your input file, the columns which need to be converted to a factor, and the column used for calculating the means
# the latter can also be a vector with multiple column names
headers <- c("chrom","pos","ctxt", "cov", "pcentmeth", "sample", "group", "organ", "sex")
toFactor <- c("sample", "group", "organ", "sex")
colsForPasteMean <- c("group")


############################################################
### from here on, there is no need to adjust anything
############################################################


# choose which model to run 
if (useDfCor == "dfCor") {
  usedModel <- dfCorrectedModel
  outfileEnding <- ".mod.dfCor"
} else {
  usedModel <- basicModel
  outfileEnding <- ".mod"
}
outfileName <- paste0(myPath, outfileEnding)

# run the models
f.print.message(paste0("processing ", myPath, "\n"))
meth <- read.table(myPath, header = FALSE, sep = "\t", col.names = headers, stringsAsFactors = FALSE)
if (max(meth$pcentmeth) <= 1) {
  f.print.message("WARNING: detected a maximum of at least 1, multiplying with 100 to get percentages\n")
  meth$pcentmeth <- 100*meth$pcentmeth
}
for (tf in toFactor) {meth[[tf]] <- as.factor(meth[[tf]])}
meth$chrPos <- paste0(meth$chrom, '_', meth$pos)
f.print.message(paste0("splitting ", myPath, "\n"))
indices <- split(1:nrow(meth), meth$chrPos)
f.print.message(paste0("testing ", length(indices), " models.\n"))
all.tests <- mclapply(indices, function(x) usedModel(meth[x,]), mc.cores = numberOfCores)
out <- do.call("rbind", all.tests)
colnames(out) <- modelOutputNames
f.print.message(paste0("writing ", outfileName, "\n"))
write.table(out, outfileName, sep = '\t', quote = FALSE)

# add the group means
f.print.message(paste0("adding group means ", outfileName, "\n"))
out <- read.table(outfileName, quote = '', sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
if (length(colsForPasteMean) > 1) {
  meth$forMean <- paste0("meanWithin_", apply(meth[,colsForPasteMean], 1, function(x) paste(x, collapse = '_')))
} else {
  meth$forMean <- meth[[colsForPasteMean]]
}
allGroups <- unique(meth$forMean)
temp <- aggregate(meth$pcentmeth, by = list(chromPos = meth$chrPos, grp = meth$forMean), function(x) mean(x, na.rm = TRUE), simplify = TRUE)
for (curGroup in allGroups) {
  subTemp <- subset(temp, grp == curGroup)
  curGroup <- paste0("meanWithin_", curGroup)
  out[[curGroup]] <- NA
  out[subTemp$chromPos, curGroup] <- subTemp$x
}
f.print.message(paste0("writing ", outfileName, "\n"))
write.table(out, outfileName, sep = '\t', quote = FALSE)
f.print.message("finished.\n")
