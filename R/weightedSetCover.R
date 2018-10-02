#' @importFrom dplyr bind_rows
#' @importFrom parallel mclapply
# Size constrained weighted set cover problem
weightedSetCover <- function(cur_config, geneset_ids, geneset_info, output_dir, n_threads){
  `%>%` <- magrittr::`%>%`
  # we only start with the top (multiplier * top_k) most
  # significant sets to reduce computational cost
  multiplier <- 10
  top_k <- cur_config[["top_num"]]
  max_num_set <- multiplier * top_k

  # sort by absolute of signP
  if(length(geneset_ids) > max_num_set){
    sorted_geneset_info <- geneset_info[order(sapply(geneset_info, function(gs) {
                                                     return(abs(gs))
                                                    }), decreasing=TRUE)]
    sorted_geneset_ids <- geneset_ids[names(sorted_geneset_info)]
    geneset_ids <- sorted_geneset_ids[1:max_num_set]
    geneset_info <- sorted_geneset_info[1:max_num_set]
  }

  s_hat <- 1.0
  # get all unique genes in all enriched sets
  all_genes <- unique(unlist(geneset_ids))
  remain <- s_hat*length(all_genes)

  # final results, contains a list of gene set names
  cur_res <- c()
  # current candidates with marginal gain and size
  all_set_names <- names(geneset_ids)
  mc_results <- mclapply(all_set_names, function(cur_name, cur_res, geneset_ids, geneset_info){
           cur_gain <- marginalGain(cur_name, cur_res, geneset_ids, geneset_info)
           cur_size <- length(geneset_ids[[cur_name]])
           return(data.frame(geneset.name=cur_name, gain=cur_gain, size=cur_size, stringsAsFactors=FALSE))
                        }, cur_res=cur_res, geneset_ids=geneset_ids, geneset_info=geneset_info, mc.cores=n_threads)
  candidates <- mc_results %>% bind_rows()
  top_k <- min(top_k, nrow(candidates))
  for(i in seq(top_k)){
    # if there is no candidates, return
    if(nrow(candidates) == 0) {
      covered_genes <- unique(unlist(geneset_ids[cur_res]))
      s_hat <- length(covered_genes)/length(all_genes)
      print("no more candidates, ending weighted set cover")
      return(list(sc_topsets=cur_res,sc_s_hat=s_hat))
    }
    # find the set with maximum marginal gain
    # tie breaker: for two sets with sname marginal gain, pick the one with
    # larger size
    candidates <- candidates[order(-candidates$gain, -candidates$size), ]
    # update remain
    remain <- remain - length(marginalBenefit(candidates[1,"geneset.name"], cur_res, geneset_ids))
    cur_res <- c(cur_res, candidates[1,"geneset.name"])
		if(remain == 0){
			covered_genes <- unique(unlist(geneset_ids[cur_res]))
			s_hat <- length(covered_genes)/length(all_genes)
			print("remain is 0, ending weighted set cover")
			# full coverage solution
			return(list(sc_topsets=cur_res,sc_s_hat=s_hat))
		}
    # update candidates
    # first remove the one just been selected
    candidates <- candidates[-1,]
    # recalculate gain, remove rows with gain == 0
    mc_results <- mclapply(seq(nrow(candidates)), function(row, candidates, cur_res, geneset_ids, geneset_info){
                             cur_name <- candidates[row, "geneset.name"]
                             cur_gain <- marginalGain(cur_name, cur_res, geneset_ids, geneset_info)
                             if(cur_gain != 0) {
                               candidates[candidates$geneset.name == cur_name, "gain"] <- cur_gain
                               tmp_candidate <- candidates[candidates$geneset.name == cur_name,]
                               return(tmp_candidate)
                             }
                        }, candidates=candidates, cur_res=cur_res, geneset_ids=geneset_ids, geneset_info=geneset_info, mc.cores=n_threads)

    new_candidates <- mc_results %>% bind_rows()
    candidates <- new_candidates
  }
  # not fully covered, compute the current coverage and return
  covered_genes <- unique(unlist(geneset_ids[cur_res]))
  s_hat <- length(covered_genes)/length(all_genes)
  print("ending weighted set cover")
  return(list(sc_topsets=cur_res,sc_s_hat=s_hat))
}

topSCGeneSets <- function(cur_config, geneset_ids, geneset_info, output_dir, n_threads){
#topSCGeneSets <- function(sc.data, top_k, n_threads){
  sc_res <- weightedSetCover(cur_config, geneset_ids, geneset_info,  output_dir, n_threads)
  save(sc_res, file=file.path(output_dir, paste0("sc_results.RData")))

  ret <- topGeneSets(cur_config, geneset_ids[sc_res$sc_topsets], geneset_info[sc_res$sc_topsets], output_dir, prefix="sc")
  return(list(sctopsets_s_hat=sc_res$sc_s_hat,
              sctopsets_num=ret$topset.num,
              sctopsets_png=ret$topset.png,
              sctopsets_pdf=ret$topset.pdf))
}

# return a list of genes from all.genes that has not been
# covered so far
# cur.set.name: name of the candidate gene set
# cur.res: vector of names of gene sets in current result
marginalBenefit <- function(cur.set.name, cur.res, geneset_ids){
  all.genes <- unique(unlist(geneset_ids))
  cur.genes <- geneset_ids[[cur.set.name]]
  if(length(cur.res) == 0){
    not.covered.genes <- cur.genes
  } else{
    covered.genes <- unique(unlist(geneset_ids[cur.res]))
    not.covered.genes <- setdiff(cur.genes, covered.genes)
  }
  return(not.covered.genes)
}


marginalGain <- function(cur.set.name, cur.res, geneset_ids, geneset_info){
  cur_abs_logp <- abs(geneset_info[[cur.set.name]])
  if(cur_abs_logp == 0) {
     stop("score cannot be zero")
  }
  cur.cost <- 1.0/cur_abs_logp
  cur.mben <- marginalBenefit(cur.set.name, cur.res, geneset_ids)
 # print(cur.cost)
 # print(cur.mben)
  return(length(cur.mben)/cur.cost)
}
