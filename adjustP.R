#!/usr/bin/env Rscript
library("stats")
library("qvalue")

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
myPath <- as.character(myarg[argPos+1])

#source("SLIMcodes.R")

cat("===", format(Sys.time(), "%Y %b %d %X"), paste("=== processing ", myPath, "\n", sep=''))
meth <- read.table(myPath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat("===", format(Sys.time(), "%Y %b %d %X"), paste("=== adjusting ", myPath, "\n", sep=''))
toAdjust <- grep("^P_", colnames(meth), value = TRUE)
for (tA in toAdjust) {
  meth[[tA]] <- as.numeric(meth[[tA]])
  tAnew <- gsub("^P", "FDR", tA)
  meth[[tAnew]] <- p.adjust(meth[[tA]], "BH")
  tAnew <- gsub("^P", "Q", tA)
  meth[[tAnew]] <- qvalue(meth[[tA]])$qvalues
  #tAnew <- gsub("^P", "QSLIM", tA)
  #slimObj <- SLIMfunc(meth[[tA]])
  #meth[[tAnew]]<- QValuesfun(meth[[tA]], slimObj$pi0_Est)
}
cat("===", format(Sys.time(), "%Y %b %d %X"), paste("=== writing ", paste0(myPath, ".adj"), "\n", sep=''))
write.table(meth, paste0(myPath, ".adj"), sep = '\t', quote = FALSE)
cat("===", format(Sys.time(), "%Y %b %d %X"), "=== finished.\n")

