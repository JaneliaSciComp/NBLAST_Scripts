if (!require("nat",character.only = TRUE)) {
  install.packages("nat", repos="http://cran.rstudio.com/")
}
if (!require("nat.nblast",character.only = TRUE)) {
  install.packages("nat.nblast", repos="http://cran.rstudio.com/")
}

library(nat)
library(nat.nblast)

#set a target directory
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) stop("At least two argument must be supplied.", call.=FALSE)
indir = args[1]
dstpath  = args[2]
indir <- gsub("\\\\", "/", indir)
dstpath <- gsub("\\\\", "/", dstpath)

cat(paste("INDIR: ", indir))
cat(paste("DSTPATH: ", dstpath))

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

label <- tools::file_path_sans_ext(basename(dstpath))
label <- paste0(label, "_")

#for swc files
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

nl <- c(nlswc, nlnrrd)

score_matrix = nblast_allbyall(nl, normalisation="mean")
saveRDS(score_matrix, file=dstpath)
print(paste("Score matrix: ", dstpath))
print(paste("SWC dir: ", file.path(dirname(dstpath), paste0(label,"swc"))))
write.neurons(nl, dir=file.path(dirname(dstpath), paste0(label,"swc")), files=names(nl), format="swc", Force=T)
sprintf("NBLAST SWCs: %i", length(nl))
write.table(format(score_matrix, digits=10), "result.txt", sep=",", quote=F, col.names=T, row.names=T)