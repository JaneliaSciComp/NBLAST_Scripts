
if (!require("nat",character.only = TRUE)) {
  install.packages("nat", repos="http://cran.rstudio.com/")
}
if (!require("nat.nblast",character.only = TRUE)) {
  install.packages("nat.nblast", repos="http://cran.rstudio.com/")
}
if (!require("foreach",character.only = TRUE)) {
  install.packages("foreach", repos="http://cran.rstudio.com/")
}
if (!require("doParallel",character.only = TRUE)) {
  install.packages("doParallel", repos="http://cran.rstudio.com/")
}


library(nat.nblast)
library(nat)
library(foreach)
library(parallel)
library(doParallel)


#set a target directory
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) stop("At least two argument must be supplied.", call.=FALSE)

nlibpath = args[1]
nlibs = strsplit(nlibpath, ",")
nlibs = nlibs[[1]]

dstpath = args[2]

cat(paste("DATABASE: ", nlibpath))
cat(paste("DSTPATH: ", dstpath))

cat("Initializing threads...\n")
th = parallel::detectCores()-1
if (th < 1) th = 1
cl <- parallel::makeCluster(th, outfile="")

tryCatch({
  
  registerDoParallel(cl)
  cat("Loading NBLAST library into each thread...\n")
  clusterCall(cl, function() suppressMessages(library(nat.nblast)))
  
  sprintf("thread num: %i", th)
  cat("Running NBLAST...\n")
  
  nl <- nlibs[[1]]
  cat(paste(nl, "\n"))
  dp <- read.neuronlistfh(nl, localdir=dirname(nl))
  
  a <- seq(1, length(dp), by=length(dp)%/%th)
  b <- seq(length(dp)%/%th, length(dp), by=length(dp)%/%th)
  if (length(dp)%%th > 0) {
    b <- c(b, length(dp))
  }
  
  cat("calculating a score matrix...\n")
  scores <- foreach(aa = a, bb = b) %dopar% {
    if (aa == 1) {
      nblast(dp[aa:bb], dp, normalised=T, UseAlpha=T, .progress='text')
    } else {
      nblast(dp[aa:bb], dp, normalised=T, UseAlpha=T)
    }
  }
  
  smat <- scores[[1]]
  if (length(scores) > 2) {
    for (i in 2:length(scores)) {
      smat <- cbind(smat, scores[[i]])
    }
  }
  
  cat("calculating mean scores...\n")
  smat <- sub_score_mat(scoremat=smat, normalisation="mean")

}, 
error = function(e) {
  print(e)
  Sys.sleep(5)
  stopCluster(cl)
})

stopCluster(cl)

cat("saving a score matrix\n")
saveRDS(smat, file=dstpath)
print(paste("Score matrix: ", dstpath))

cat("Done\n")