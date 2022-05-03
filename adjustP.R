#!/usr/bin/env Rscript
library("stats")

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
myPath <- as.character(myarg[argPos+1])
adjType <- as.character(myarg[argPos+2])

if (adjType == "Q") {
  source("SLIMcodes.R")
} else if (adjType == "QSLIM") {
  library("qvalue")
} else {
  adjType <- "FDR"
}

f.print.message <- function(x) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste0("=== ", x,"\n")) }

f.print.message(paste0("processing ", myPath, " (", adjType, ")"))
meth <- read.table(myPath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
f.print.message(paste0("adjusting ", myPath))
toAdjust <- grep("^P_", colnames(meth), value = TRUE)
for (tA in toAdjust) {
  meth[[tA]] <- as.numeric(meth[[tA]])
  if (adjType == "Q") {
    tAnew <- gsub("^P", "Q", tA)
    meth[[tAnew]] <- qvalue(meth[[tA]])$qvalues
  } else if (adjType == "QSLIM") {
    tAnew <- gsub("^P", "QSLIM", tA)
    slimObj <- SLIMfunc(meth[[tA]])
    meth[[tAnew]]<- QValuesfun(meth[[tA]], slimObj$pi0_Est)
  } else {
    tAnew <- gsub("^P", "FDR", tA)
    meth[[tAnew]] <- p.adjust(meth[[tA]], "BH")
  }
}
f.print.message(paste0("writing ", paste0(myPath, ".adj")))
write.table(meth, paste0(myPath, ".adj"), sep = '\t', quote = FALSE)
f.print.message("finished.")


