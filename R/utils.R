# convert list of gene names to gmt file
# each list component has a gene set name
# the value of each list item is an array of gene names/ids
list2gmt <- function(rlist, out_file) {
  m <- do.call(rbind, lapply(seq_along(rlist), function(i) paste(c(names(rlist)[[i]],"", rlist[[i]]), collapse="\t")))
  write.table(m,file=out_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# for now, for my own use
extractScore <- function(rlist, out_file) {
  m <- do.call(rbind, lapply(seq_along(rlist), function(i) paste(c(names(rlist)[[i]], rlist[[i]][["signP"]]), collapse="\t")))
  write.table(m,file=out_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
}


gmt2list <- function(gmt_file){
  if (!file.exists(gmt_file)) {
    stop("There is no such gmt file.")
  }
  x <- scan(gmt_file, what="", sep="\n", quiet=TRUE)
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  gmt_list <- lapply(y, `[`, c(-1,-2))
}

# for cytoscape, set node size
computeNodeDim <- function(maxVal, minVal, val){
  maxCyNodeSize<-200;
  minCyNodeSize<-50;
  return((val-minVal)/(maxVal-minVal)*(maxCyNodeSize-minCyNodeSize)+minCyNodeSize)
}


# set color: min: blue, middle: white, max: red
computeNodeColor <- function(maxVal, minVal, val){
  # if all positive, map to [0.5, 1]
  if(minVal > 0){
    maxColor <- 1;
    minColor <- 0.5;
  } else if(maxVal < 0) {  # all negative
    maxColor <- 0.5;
    minColor <- 0;
  } else {
    larger = max(maxVal, abs(minVal));
    maxVal <- larger;
    minVal <- -1.0*larger;
    maxColor <- 1;
    minColor <- 0;
  }
  return((val-minVal)/(maxVal-minVal)*(maxColor-minColor)+minColor)
}
