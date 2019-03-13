
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) stop("At least two argument must be supplied.", call.=FALSE)
indir = args[1]

nlibpath  = args[2]
nlibs = strsplit(nlibpath, ",")
nlibs = nlibs[[1]]

if (length(args) >= 3) {
  outfname = args[3]
  if (outfname == "none") outfname = ""
} else {
  outfname = paste(basename(imagefile), ".nblust", sep="")
}

if (length(args) >= 4) {
	outputdir = args[4]
} else {
	outputdir = dirname(imagefile)
}
if (!dir.exists(outputdir)) {
    dir.create(outputdir, FALSE)
}

if (length(args) >= 5) {
	resultnum = strtoi(args[5])
	if (is.na(resultnum)) resultnum = 10
} else {
  resultnum = 10
}

if (length(args) >= 6) {
  dbnames = args[6]
} else {
  dbnames = paste(basename(nlibs), collapse=",")
}

if (length(args) >= 7) {
  normalization = args[7]
} else {
  normalization = "mean"
}

if (length(args) >= 8) {
  rsmp = as.double(args[8])
  if (is.na(rsmp)) rsmp = 3.0
} else {
  rsmp = 3.0
}

if (length(args) >= 9) {
  kval = strtoi(args[9])
  if (is.na(kval)) kval = 3
} else {
  kval = 3
}

if (length(args) >= 10) {
  thnum = strtoi(args[10])
  if (is.na(thnum)) thnum = 0
} else {
  thnum = 0
}

cat("Loading NAT...")

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

#read input neurons
#swc files
swclist <- list.files(indir, pattern='swc$', full.names=F)
swclist <- file.path(indir, swclist)
sprintf("SWCs: %i", length(swclist))
nlswc <- neuronlist()
if (length(swclist) > 0) {
  neu = read.neurons(swclist)
  nlswc = dotprops(neu, resample=rsmp, OmitFailures=T, k=kval)
  fnames <- names(nlswc)
  nnames <- list()
  for (i in 1:length(fnames)) {
    nnames[[i]] <- tools::file_path_sans_ext(fnames[i])
  }
  names(nlswc) <- nnames
  sprintf("swc neurons: %i", length(nlswc))
}

#nrrd files (skeletonized images)
nlnrrd <- neuronlist()
nrrdlist <- list.files(indir, pattern='nrrd$', full.names=F)
sprintf("NRRDs: %i", length(nrrdlist))
if (length(nrrdlist) > 0) {
  for (i in 1:length(nrrdlist)) {
    tryCatch({
      img <- read.im3d(file.path(indir, nrrdlist[i]))
      dp <- dotprops(img, OmitFailures=T, k=kval)
      nlnrrd[[i]] <- dp
      print(i)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  nnames2 <- list()
  for (i in 1:length(nrrdlist)) {
    nnames2[[i]] <- tools::file_path_sans_ext(nrrdlist[i])
  }
  names(nlnrrd) <- nnames2
  sprintf("nrrd neurons: %i", length(nlnrrd))
}

innl <- c(nlswc, nlnrrd)


cat("Initializing threads...\n")
thmax = parallel::detectCores()-1
if ( (thnum < 1) || (thnum > thmax) ) thnum = thmax
cl <- parallel::makeCluster(thnum, outfile="")


tryCatch({
  
  registerDoParallel(cl)
  cat("Loading NBLAST library into each thread...\n")
  clusterCall(cl, function() suppressMessages(library(nat.nblast)))
  
  sprintf("thread num: %i", thnum)
  cat("Running NBLAST...\n")
  
  for (nid in 1:length(innl)) {
    qdp <- innl[[nid]]
    
    allres <- neuronlist();
    allscr <- numeric();
  
    for (i in 1:length(nlibs)) {
      nl <- nlibs[i]
      cat(paste(nl, "\n"))
      dp <- read.neuronlistfh(nl, localdir=dirname(nl))
      
      a <- seq(1, length(dp), by=length(dp)%/%thnum)
      b <- seq(length(dp)%/%thnum, length(dp), by=length(dp)%/%thnum)
      if (length(dp)%%thnum > 0) {
        b <- c(b, length(dp))
      }
      
      cat("calculating forward scores...\n")
      fwdscores <- foreach(aa = a, bb = b) %dopar% {
        if (aa == 1) {
          nblast(qdp, dp[aa:bb], normalised=T, UseAlpha=T, .progress='text')
        } else {
          nblast(qdp, dp[aa:bb], normalised=T, UseAlpha=T)
        }
      }
      scnames <- list()
      for (j in 1:length(fwdscores)) scnames <- c(scnames, names(fwdscores[[j]]))
      fwdscores <- unlist(fwdscores)
      names(fwdscores) <- scnames
      
      if (normalization == "mean") {
        cat("calculating reverse scores...\n")
        revscores <- foreach(aa = a, bb = b) %dopar% {
          if (aa == 1) {
            nblast(dp[aa:bb], qdp, normalised=T, UseAlpha=T, .progress='text')
          } else {
            nblast(dp[aa:bb], qdp, normalised=T, UseAlpha=T)
          }
        }
        scnames <- list()
        for (j in 1:length(revscores)) scnames <- c(scnames, names(revscores[[j]]))
        revscores <- unlist(revscores)
        names(revscores) <- scnames
        
        scores <- (fwdscores + revscores) / 2
      } else {
        scores <- fwdscores
      }
      
      cat("sorting scores...\n")
      scores <- sort(scores, dec=T)
      
      if (length(scores) <= resultnum) {
        #results <- as.neuronlist(dp[names(scores)])
        slist <- scores
      } else {
        #results <- as.neuronlist(dp[names(scores)[1:resultnum]])
        slist <- scores[1:resultnum]
      }
      
      cat("setting names...\n")
      #for (j in 1:length(results)) {
      #  names(results)[j] <- paste(names(results[j]), as.character(i-1), sep=",")
      #}
      for (j in 1:length(slist)) {
        names(slist)[j] <- paste(names(slist[j]), as.character(i-1), sep=",")
      }
      
      cat("combining lists...\n")
      #allres <- c(allres, results)
      allscr <- c(allscr, slist)
      
      rm(dp)
      gc()
    }
    
    allscr = sort(allscr, dec=T)
    if (length(allscr) <= resultnum) {
      #results = allres[names(allscr)]
      slist = allscr
    } else {
      #results = allres[names(allscr)[1:resultnum]]
      slist = allscr[1:resultnum]
    }
    
    cat("Writing results...\n")
    outfname <- names(innl)[[nid]]
    #swczipname <- paste(outfname, ".zip", sep="")
    rlistname  <- paste(outfname, ".txt", sep="")
    #zprojname  <- paste(outfname, ".png", sep="")
    
    f = file(file.path(outputdir,rlistname))
    writeLines(c(dbnames, nlibpath), con=f)
    write.table(format(slist, digits=8, scientific = FALSE), append=T, file.path(outputdir,rlistname), sep=",", quote=F, col.names=F, row.names=T)
    #n = names(results)
    #n = gsub(",", " db=", n)
    #write.neurons(results, dir=file.path(outputdir,swczipname), files=n, format='swc', Force=T)
    
    #if (!is.null(img)) {
    #  zproj = projection(img, projfun=max)
    #  size = dim(zproj)
    #  png(file.path(outputdir,zprojname), width=size[1], height=size[2])
    #  par(plt=c(0,1,0,1))
    #  image(zproj, col = grey(seq(0, 1, length = 256)))
    #  dev.off()
    #} else {
    #  png(file.path(outputdir,zprojname), width=256, height=256, bg='black')
    #  plot.new()
    #  dev.off()
    #}
  }
  
}, 
error = function(e) {
  print(e)
  Sys.sleep(5)
  stopCluster(cl)
})

stopCluster(cl)

cat("Done\n")
