#' affinityPropagation
#'
#' @param apdata:  a list of input data
#'       $ output.dir
#'       $ genesetInfo
#'       $ genesetIds
#' @param sim method for compute similarity between gene set:
#'  default is 'Jaccard'
#' @importFrom apcluster apcluster
#' @importFrom apcluster aggExCluster
#' @importFrom jsonlite unbox toJSON
#' @import proxy
affinityPropagation <- function(ap.data, sim="Jaccard"){
  output.dir <- ap.data$output.dir
  genesetIds <- ap.data$genesetIds
  genesetInfo <- ap.data$genesetInfo
  md5val <- ap.data$md5val

  # compute the similiarity and input preference vector
  if(sim == "Jaccard"){
    ret <- jaccardSim(genesetInfo, genesetIds)
  } else {
    stop("similarity algorithm not supported; choose one of the following:
         Jaccard;")
  }
  sim.mat <- ret$sim.mat
  ip.vec <- ret$ip.vec

  ap.result <- apcluster(sim.mat, p=ip.vec)

  # create json file for cytoscape
  json.data <- list()
  vis.data <- list()
  clusters <- ap.result@clusters
  exemplars <- ap.result@exemplars

  # run affinity propagation among exemplars if there are more
  # than 10 clusters
  clusters2 <- NULL
  exemplars2 <- NULL
  ap.result2 <- NULL
  sim.mat2 <- NULL
  min.num.clusters <- 10
  genesetIds2 <- genesetIds[names(exemplars)]
  genesetInfo2 <- genesetInfo[names(exemplars)]
  if(length(clusters) > min.num.clusters){
    if(sim == "Jaccard"){
      ret2 <- jaccardSim(genesetInfo2, genesetIds2)
    }
    sim.mat2 <- ret2$sim.mat
    ap.result2 <- apcluster(sim.mat2)
    clusters2 <- ap.result2@clusters
    exemplars2 <- ap.result2@exemplars
  }

  cnt <- 1
  node.df <- data.frame(name=character(), platform=character(),
                        size=numeric(), color=numeric())
  # get the max, min value of leNum for all nodes
  all.leNum <- sapply(genesetInfo, function(x) { return(x$leNum) })
  all.NES <- sapply(genesetInfo, function(x) { return(x$NES) })
  all.platforms <- sapply(genesetInfo, function(x) { return(x$platform) })
  all.platforms <- unique(all.platforms)
  maxLeNum <- max(all.leNum)
  minLeNum <- min(all.leNum)
  maxNES <- max(all.NES)
  minNES <- min(all.NES)
  # nodes
  for(i in seq_along(clusters)){
    nodes <- names(clusters[[i]])
    cur.cnt <- 1
    for(j in seq_along(nodes)){
      dat <- list()
      dat[["id"]] <- unbox(nodes[j])
      if(nodes[j] %in% names(exemplars)){
        dat[["is_exemplar"]]  <- unbox(1)
      } else{
        dat[["is_exemplar"]]  <- unbox(0)
      }
      if(!is.null(exemplars2)){
        if(nodes[j] %in% names(exemplars2)){
          dat[["is_exemplar"]]  <- unbox(2)
        }
      }
      dat[["signP"]] <- unbox(genesetInfo[[nodes[j]]][["signP"]])
      dat[["NES"]] <- unbox(genesetInfo[[nodes[j]]][["NES"]])
      dat[["nSize"]] <- unbox(computeNodeDim(maxLeNum, minLeNum, genesetInfo[[nodes[j]]][["leNum"]]))
      dat[["color"]] <- unbox(computeNodeColor(maxNES, minNES, genesetInfo[[nodes[j]]][["NES"]]))
      dat[["platform"]] <- unbox(which(all.platforms == genesetInfo[[nodes[j]]][["platform"]]))
      vis.data[[cnt]] <- list(data=dat)
      cnt <- cnt+1
      node.df <- rbind(node.df, list(name=nodes[j], platform=dat[["platform"]],
                                     size=dat[["nSize"]], color=dat[["color"]]), stringsAsFactors=FALSE)
      cur.cnt <- cur.cnt+1
    }
  }
  node_file <- file.path(output.dir, paste0("ap_", md5val,"_nodelist.txt"))
  write.table(node.df, node_file, col.names = TRUE, sep='\t',quote=FALSE, row.name=FALSE)

  # edges
  edge.cnt <- 1
  edge.df <- data.frame(node1=character(), node2=character())
  for(j in seq_along(exemplars)){
    cur.exemplar <- names(exemplars[j])
    cur.cluster <- names(clusters[[j]])
    for(k in seq_along(cur.cluster)){
      if(cur.exemplar != cur.cluster[k]){
        # edge level:
        # level 1:  between an exemplar and a node
        # level 2: between exemplar of exemplar and original exemplar
        # for large networks, we may only show level 2 edges
        vis.data[[cnt]] <- list(data=list(id=unbox(paste0("e",edge.cnt)), level=unbox(1),
                                          source=unbox(cur.exemplar), target=unbox(cur.cluster[k])))
        edge.df <- rbind(edge.df, list(node1=cur.exemplar, node2=cur.cluster[k]), stringsAsFactors=FALSE)
        cnt <- cnt+1
        edge.cnt <- edge.cnt+1
      }
    }
  }
  if(!is.null(exemplars2)){
    for(j in seq_along(exemplars2)){
      cur.exemplar <- names(exemplars2[j])
      cur.cluster <- names(clusters2[[j]])
      for(k in seq_along(cur.cluster)){
        if(cur.exemplar != cur.cluster[k]){
          # edge level:
          # level 1:  between an exemplar and a node
          # level 2: between exemplar of exemplar and original exemplar
          # for large networks, we may only show level 2 edges
          vis.data[[cnt]] <- list(data=list(id=unbox(paste0("e",edge.cnt)), level=unbox(2),
                                            source=unbox(cur.exemplar), target=unbox(cur.cluster[k])))
          edge.df <- rbind(edge.df, list(node1=cur.exemplar, node2=cur.cluster[k]), stringsAsFactors=FALSE)
          cnt <- cnt+1
          edge.cnt <- edge.cnt+1
        }
      }
    }
  }
  json.data[["vis"]] <- vis.data
  json_file <- file.path(output.dir, paste0("ap_",md5val,".json"))
  write(toJSON(json.data), json_file)
  edge_file <- file.path(output.dir, paste0("ap_", md5val,"_edgelist.txt"))
  write.table(edge.df, edge_file, col.names = FALSE, sep='\t',quote=FALSE, row.name=FALSE)
  # create a zip file for download
  zip_file <- file.path(output.dir, paste0("ap_",md5val,".zip"))
  zip(zip_file, c(node_file, edge_file, json_file), flags="-j")
  return(list(json=basename(json_file), zip=basename(zip_file)))
}


jaccardSim <- function(genesetInfo, genesetIds){
  # first find out the union of sets, sorted
  all.genes <- sort(unique(unlist(genesetIds)))
  overlap.mat <- sapply(genesetIds, function(x) {as.integer(all.genes %in% x)})
  # proxy::simil
  sim.mat <- as.matrix(simil(overlap.mat, by_rows=FALSE, method="Jaccard"))
  sim.mat[is.na(sim.mat)] <- 1
  # if there is no overlap, set the similarity to -Inf
  sim.mat[sim.mat == 0] <- -Inf
  # check sim.mat to see if it is identical for each pair
  if(max(sim.mat)==min(sim.mat)){
    # this will generate error, so randomy add some noise to off diagonal elements
    mat.siz <- dim(sim.mat)[1]
    rand.m <- matrix(rnorm(mat.siz*mat.siz,0,0.01),mat.siz)
    # make it symmetric
    rand.m[lower.tri(rand.m)] = t(rand.m)[lower.tri(rand.m)]
    sim.mat <- sim.mat + rand.m
    # make diagonal all 1
    diag(sim.mat) <- 1
  }
  # set the input preference (IP) for each geneset
  # give higher IP to geneset with larger -logP (remove sign)
  # IP <- maxScore for geneset with largest -logP value
  # IP <- minScore for geneset with smallest -logP value
  # other genesets will have linearly interpolated IP value
  all.minus.logP <- sapply(genesetInfo, function(x){
    return(abs(x$signP))
  })
  max.minus.logP <- max(all.minus.logP)
  min.minus.logP <- min(all.minus.logP)
  # for Jaccard
  #minScore <- min(sim.mat)
  minScore <- 0
  #maxScore <- max(sim.mat)
  tmp.sim.mat <- sim.mat
  tmp.sim.mat[!is.finite(tmp.sim.mat)] <- NA
  # get the median excluding -Inf
  maxScore <- median(tmp.sim.mat, na.rm=TRUE)
  if(abs(max.minus.logP-min.minus.logP) < .Machine$double.eps^0.5){
    ip.vec <- NA
  } else{
    ip.vec <- sapply(genesetInfo, function(x){
      return(minScore+(maxScore-minScore)*(abs(x$signP)-min.minus.logP)/(max.minus.logP-min.minus.logP))
    })
  }
  return(list(sim.mat=sim.mat, ip.vec=ip.vec))
}
