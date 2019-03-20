
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



seglength=function(ThisSeg, sum=TRUE){
  #ThisSeg is an array of x,y and z data points
  #In order to calculate the length
  #Need to find dx,dy,dz
  #Then find sqrt(dx^2+...)
  #Then sum over the path
  if(nrow(ThisSeg)==1) return(0)
  ds=diff(ThisSeg)
  edgelengths=sqrt(rowSums(ds*ds))
  if(sum) sum(edgelengths) else unname(edgelengths)
}

resample2.neuron<-function(x, stepsize, ...) {
  # extract original vertex array before resampling
  cols=c("X","Y","Z")
  if(!is.null(x$d$W)) cols=c(cols, 'W')
  # if(!is.null(x$d$Label)) cols=c(cols, 'Label')
  d=data.matrix(x$d[, cols, drop=FALSE])
  # if(!is.null(d$Label)) d$Label=as.integer(d$Label)
  if(any(is.na(d[,1:3])))
    stop("Unable to resample neurons with NA points")
  
  # fetch all segments and process each segment in turn
  sl=as.seglist(x, all = T, flatten = T)
  npoints=nrow(d)
  dl=list(d)
  sl2=list()
  count=0
  for (i in seq_along(sl)){
    s=sl[[i]]
    # interpolate this segment
    dold=d[s, , drop=FALSE]
    dnew=resample_segment2(dold, stepsize=stepsize, ...)
    if(is.null(dnew)) next
    dl[[length(dl)+1]]=dnew
    # if we've got here, we need to do something
    # add new points to the end of the swc block
    # and give them sequential point numbers
    newids=seq.int(from = npoints+1, length.out = nrow(dnew))
    npoints=npoints+nrow(dnew)
    # replace internal ids in segment so that proximal point is connected to head
    # and distal point is connected to tail
    sl[[i]]=c(s[1], newids, s[length(s)])
    
    count = count + 1
    sl2[[count]] = sl[[i]]
  }
  sl = sl2
  if (length(sl) == 0) {
    return(NULL)
  }
  
  d=do.call(rbind, dl)
  d=as.data.frame(d)
  rownames(d)=NULL
  # let's deal with the label column which was dropped - assume that always the
  # same within a segment
  head_idxs=sapply(sl, "[", 1)
  seglabels=x$d$Label[head_idxs]
  
  # in order to avoid re-ordering the segments when as.neuron.ngraph is called
  # we can renumber the raw indices in the seglist (and therefore the vertices)
  # in a strictly ascending sequence based on the seglist
  # it is *much* more efficient to compute this on a single vector rather than
  # separately on each segment in the seglist. However this does involve some
  # gymnastics 
  usl=unlist(sl)
  old_ids=unique(usl)
  # reorder vertex information to match this
  d=d[old_ids,]
  
  node_per_seg=sapply(sl, length)
  df=data.frame(id=usl, seg=rep(seq_along(sl), node_per_seg))
  df$newid=match(df$id, old_ids)
  sl=split(df$newid, df$seg)
  labels_by_seg=rep(seglabels, node_per_seg)
  # but there will be some duplicated ids (branch points) that we must remove
  d$Label=labels_by_seg[!duplicated(df$newid)]
  swc=seglist2swc(sl, d)
  as.neuron(swc)
}

# Interpolate ordered 3D points (optionally w diameter)
# NB returns NULL if unchanged (when too short or <=2 points) 
# and only returns _internal_ points, omitting the head and tail of a segment
#' @importFrom stats approx
resample_segment2<-function(d, stepsize, ...) {
  # we must have at least 2 points to resample
  if(nrow(d) < 2) return(NULL)
  
  dxyz=xyzmatrix(d)
  # we should only resample if the segment is longer than the new stepsize
  l=seglength(dxyz)
  if(l<=stepsize) return(NULL)
  
  # figure out linear position of new internal points
  internalPoints=seq(stepsize, l, by=stepsize)
  nInternalPoints=length(internalPoints)
  # if the last internal point is actually in exactly the same place 
  # as the endpoint then discard it
  if(internalPoints[nInternalPoints]==l) {
    internalPoints=internalPoints[-length(internalPoints)]
    nInternalPoints=length(internalPoints)
  }
  
  # find cumulative length stopping at each original point on the segment
  diffs=diff(dxyz)
  cumlength=c(0,cumsum(sqrt(rowSums(diffs*diffs))))
  
  # find 3D position of new internal points
  # using linear approximation on existing segments
  # make an emty object for results
  # will have same type (matrix/data.frame as input)
  dnew=matrix(nrow=nInternalPoints, ncol=ncol(d))
  colnames(dnew)=colnames(d)
  if(is.data.frame(d)){
    dnew=as.data.frame(dnew)
  }
  for(n in seq.int(ncol(dnew))) {
    dnew[,n] <- if(!all(is.finite(d[,n]))) {
      rep(NA, nInternalPoints)
    } else {
      approx(cumlength, d[,n], internalPoints, 
             method = ifelse(is.double(d[,n]), "linear", "constant"))$y
    }
  }
  dnew
}



#read input neurons
#swc files
swclistfn <- list.files(indir, pattern='swc$', full.names=F)
swclist <- file.path(indir, swclistfn)
sprintf("SWCs: %i", length(swclist))
nlswc <- neuronlist()
if (length(swclist) > 0) {
  for (i in 1:length(swclist)) {
    tryCatch({
      n = read.neuron(swclist[[i]])
      rn = NULL
      step = rsmp
      repeat{
        rn = resample2.neuron(n, stepsize = step)
        if (!is.null(rn) || step <= 0.5) break
        step = step / 2.0
      }
      if (!is.null(rn)) {
        dp = dotprops(rn, OmitFailures=T, k=kval)
        nlswc[[i]] <- dp
      } else {
        cat("ERROR : the stepsize is too large.\n")
        cat("Path :", swclist[[i]], "\n")
      }
    }, error=function(e){
      cat("ERROR :", conditionMessage(e), "\n")
      cat("Path :", swclist[[i]], "\n")
    })
  }
  nnames <- list()
  for (i in 1:length(swclistfn)) {
    nnames[[i]] <- tools::file_path_sans_ext(swclistfn[[i]])
  }
  names(nlswc) <- nnames
  nlswc[sapply(nlswc, is.null)] <- NULL
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
    }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
  }
  nnames2 <- list()
  for (i in 1:length(nrrdlist)) {
    nnames2[[i]] <- tools::file_path_sans_ext(nrrdlist[[i]])
  }
  names(nlnrrd) <- nnames2
  nlnrrd[sapply(nlnrrd, is.null)] <- NULL
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
    
    allfpath <- list();
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
        slist <- scores
      } else {
        slist <- scores[1:resultnum]
      }
      
      fpaths <- list()
      libmipdir <- file.path(dirname(nl), "swc_prev")
      libswcdir <- file.path(dirname(nl), "swc")
      for (j in 1:length(slist)) {
        mippath <- file.path(libmipdir, paste0(names(slist[j]), ".png"))
        swcpath <- file.path(libswcdir, paste0(names(slist[j]), ".swc"))
        tmp <- list(mippath, swcpath)
        names(tmp) <- c("mip", "swc")
        fpaths[[j]] <- tmp
      }
      
      cat("setting names...\n")
      for (j in 1:length(slist)) {
        names(fpaths)[j] <- paste(names(slist[j]), as.character(i-1), sep=",")
        names(slist)[j] <- paste(names(slist[j]), as.character(i-1), sep=",")
      }
      
      cat("combining lists...\n")
      allfpath <- c(allfpath, fpaths)
      allscr <- c(allscr, slist)
      
      rm(dp)
      gc()
    }
    
    allscr = sort(allscr, dec=T)
    if (length(allscr) <= resultnum) {
      datapaths = allfpath[names(allscr)]
      slist = allscr
    } else {
      datapaths = allfpath[names(allscr)[1:resultnum]]
      slist = allscr[1:resultnum]
    }
    
    cat("Writing results...\n")
    outfname <- names(innl)[[nid]]
    rlistname  <- paste(outfname, ".txt", sep="")

    f = file(file.path(outputdir,rlistname))
    writeLines(c(dbnames, nlibpath), con=f)
    write.table(format(slist, digits=8, scientific = FALSE), append=T, file.path(outputdir,rlistname), sep=",", quote=F, col.names=F, row.names=T)
    
    
    result_dir <- file.path(outputdir, outfname)
    dir.create(result_dir, recursive = TRUE)
    strw <- nchar( toString(format(resultnum, scientific = F)) )
    for (i in 1:length(slist)) {
      tryCatch({
        
        rank <- formatC(i, width = strw, format = "d", flag = "0")
        intscore <- as.integer(slist[[i]]*10000)
        score <- formatC(intscore, width = 5, format = "d", flag = "0")
        header <- paste(rank, score, sep="_")
        
        srcmippath <- datapaths[[i]][["mip"]]
        dstmippath <- file.path(result_dir, paste(header, basename(srcmippath), sep="_") ) 
        if (file.exists(srcmippath)) {
          file.copy(srcmippath, dstmippath, overwrite=TRUE)
        }
        
        srcswcpath <- datapaths[[i]][["swc"]]
        dstswcpath <- file.path(result_dir, paste(header, basename(srcswcpath), sep="_"))
        if (file.exists(srcswcpath)) {
          file.copy(srcswcpath, dstswcpath, overwrite=TRUE)
        }
        
      }, error=function(e){
        cat("ERROR :", conditionMessage(e), "\n")
      })
    }
    
  }
  
}, 
error = function(e) {
  print(e)
  Sys.sleep(5)
  stopCluster(cl)
})

stopCluster(cl)

cat("Done\n")
