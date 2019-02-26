library(nat)
library(nat.nblast)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) stop("At least two argument must be supplied.", call.=FALSE)
indir = args[1]
dstdb  = args[2]
indir <- gsub("\\\\", "/", indir)
indir <- gsub("/$", "", indir)
dstdb <- gsub("\\\\", "/", dstdb)

if (!file.exists(dstdb)) stop("database not found", call.=FALSE)

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

nladd <- c(nlswc, nlnrrd)

#load dst database
dstddir <- file.path(dirname(dstdb), "data")
nlfhdst <- read.neuronlistfh(dstdb , localdir=dstddir)
dstkm <- attr(nlfhdst, "keyfilemap")

#delete duplicates
idnums <- which( names(dstkm) %in% names(nladd) )
for (i in 1:length(idnums)) {
  fname <- file.path( dstddir, dstkm[[ idnums[i] ]] )
  file.remove(fname)
  print(paste("remove dup: ", fname))
}
if (length(idnums) > 0) {
  delids <- sapply(idnums, "*", -1)
  dstkm <- dstkm[delids]
}

#add neurons into dst database
nlfhadd <- as.neuronlistfh(nladd, dbdir=dstddir)
addkm <- attr(nlfhadd, "keyfilemap")

rkm <- c(dstkm, addkm)
db=filehash::dbInit(dstddir, type='RDS')
res <- neuronlistfh(db, NULL, rkm)

#overwrite
write.neurons(nladd, dir=file.path(dirname(dstdb), "swc"), files=names(nladd), format="swc", Force=T)
write.neuronlistfh(res, file=dstdb, overwrite=TRUE)
