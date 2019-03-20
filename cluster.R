if (!require("nat",character.only = TRUE)) {
  install.packages("nat", repos="http://cran.rstudio.com/")
}
if (!require("nat.nblast",character.only = TRUE)) {
  install.packages("nat.nblast", repos="http://cran.rstudio.com/")
}
if (!require("dendroextras",character.only = TRUE)) {
  install.packages("dendroextras", repos="http://cran.rstudio.com/")
}
if (!require("randomcoloR",character.only = TRUE)) {
  install.packages("randomcoloR", repos="http://cran.rstudio.com/")
}

library(nat)
library(nat.nblast)
library(dendroextras)
library(randomcoloR)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) stop("At least four argument must be supplied.", call.=FALSE)
scmatpath = args[1]
swcd = args[2]
mipd = args[3]
od  = args[4]
scmatpath <- gsub("\\\\", "/", scmatpath)
swcd <- gsub("\\\\", "/", swcd)
mipd <- gsub("\\\\", "/", mipd)
od <- gsub("\\\\", "/", od)

if (length(args) >= 5) {
  mtd = args[5]
} else {
  mtd = "ward.D"
}

if (length(args) >= 6) {
  kval = strtoi(args[6])
  if (is.na(kval)) kval = 50
} else {
  kval = 50
}

if (length(args) >= 7) {
  hval = as.double(args[7])
  if (is.na(hval)) hval = 0.0
} else {
  hval = 0.0
}

score_matrix <- readRDS(scmatpath)
hckcs <- nhclust(scoremat=score_matrix, method=mtd, maxneurons=500000)
if (hval <= 0.0) {
  sl <- slice(hckcs, k=kval)
  dkcs <- colour_clusters(hckcs, k=kval, col=randomColor(length(unique(sl)), hue = "random", luminosity="bright"))
} else {
  sl <- slice(hckcs, h=hval)
  dkcs <- colour_clusters(hckcs, h=hval, col=randomColor(length(unique(sl)), hue = "random", luminosity="bright"))
}

print(length(unique(sl)))

newlabel <- list()
for (i in 1:length(labels(dkcs))) {
  str <- paste(labels(dkcs)[[i]], sl[[labels(dkcs)[[i]]]])
  newlabel[[i]] <- str
}
labels(dkcs) <- newlabel

svgw <- 15.0*(nrow(score_matrix)/50.0)
if (svgw < 15.0) svgw <- 15.0
svgh <- 15.0*(nrow(score_matrix)/50.0)
if (svgh < 15.0) svgh <- 15.0

pdf(file.path(od, "dendrogram.pdf"), width=svgw, height=svgh/2.0, version="1.7")
plot(dkcs, xlab="", ylab="", main="", sub="", axes=TRUE)
dev.off()

if (dir.exists(swcd) || dir.exists(mipd)) {
  clusters_dir <- file.path(od,"clusters")
  strw <- nchar( toString(format(length(unique(sl)), scientific = F)) )
  for (i in 1:length(sl)) {
    tryCatch({
      dirname <- paste("cluster_", formatC(sl[[i]], width = strw, format = "d", flag = "0"), sep="")
      cdir <- file.path( clusters_dir, dirname )
      dir.create(cdir, recursive = TRUE)
      swcname <- paste(names(sl)[[i]], ".swc", sep="")
      dstpath <- file.path(cdir, swcname)
      srcpath <- file.path(swcd, swcname)
      if (file.exists(srcpath)) {
        file.copy(srcpath, dstpath, overwrite=TRUE)
      }
      pngname <- paste(names(sl)[[i]], ".png", sep="")
      dstpath <- file.path(cdir, pngname)
      srcpath <- file.path(mipd, pngname)
      if (file.exists(srcpath)) {
        file.copy(srcpath, dstpath, overwrite=TRUE)
      }
    }, error=function(e){
      cat("ERROR :", conditionMessage(e), "\n")
    })
  }
}