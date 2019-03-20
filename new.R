library(nat)
library(nat.nblast)

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


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) stop("At least two argument must be supplied.", call.=FALSE)
indir = args[1]
dstdb  = args[2]
indir <- gsub("\\\\", "/", indir)
indir <- gsub("/$", "", indir)
dstdb <- gsub("\\\\", "/", dstdb)

if (length(args) >= 3) {
  rsmp = as.double(args[3])
  if (is.na(rsmp)) rsmp = 3.0
} else {
  rsmp = 3.0
}

if (length(args) >= 4) {
  kval = strtoi(args[4])
  if (is.na(kval)) kval = 3
} else {
  kval = 3
}

print(paste("INDIR: ", indir))
print(paste("DSTDB: ", dstdb))
sprintf("RESAMPLE: %f", rsmp)
sprintf("K_VALUE: %i", kval)

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

nlnew <- c(nlswc, nlnrrd)
sprintf("total neurons: %i", length(nlnew))

#create new database
dstddir <- file.path(dirname(dstdb), "data")
nlfhnew <- as.neuronlistfh(nlnew, dbdir=dstddir)
if (length(nlnrrd) > 0) write.neurons(nlnrrd, dir=file.path(dirname(dstdb), "swc"), files=names(nlnrrd), format="swc", Force=T)
write.neuronlistfh(nlfhnew, file=dstdb, overwrite=TRUE)
