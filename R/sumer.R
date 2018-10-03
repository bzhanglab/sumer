#' summarizing multiomics enrichment analysis results
#'
#' @param config_file: specifiy input data locations
#' @param output_dir: output directory
#' @export
sumer <- function(config_file, output_dir, n_threads=4){
  full_config <- check_config(config_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  output_index_file <- file.path(output_dir, "index.html")
  conn <- file(system.file("assets", "head.html", package="sumer"), open="r")
  lines <- readLines(conn)
  for (i in 1:length(lines)){
    cat(lines[i], file=output_index_file, append=TRUE, sep="\n")
  }
  close(conn)

  cat("\n<body>\n",file=output_index_file,append=TRUE)
  cat('<div id="top" style="margin-bottom:0;">\n',file=output_index_file,append=TRUE)
  cat('<iframe src="https://s3-us-west-2.amazonaws.com/sumer-r-app/top.html" style="border:0px;height:150px;width:100%"></iframe>\n</div>\n',file=output_index_file,append=TRUE)
  cat('<div style="margin-left:10px" class="container">\n', file=output_index_file, append=TRUE, sep='')
  cat('<ol class="breadcrumb breadcrumb-arrow" style="margin:0;">\n', file=output_index_file, append=TRUE, sep='')
  cat('<li><a href="#">', full_config$project, '</a></li>\n', file=output_index_file, append=TRUE, sep='')
  cat('</ol>\n</div>\n', file=output_index_file, append=TRUE, sep='')

  config <- full_config$data
  n_platform <- nrow(config)
  all_platform_abbr <- c()

  # plot the original data
  all_data <- data.frame(platfrom=character(), name=character(), size=integer(), score=numeric)

  # save the name of output directory for each platform
  work_dirs <- list()
  for(i in seq_len(n_platform)){
    work_dir <- uuid.gen()
    cur_output_dir <- file.path(output_dir, work_dir)
    if(!dir.exists(cur_output_dir)) {
      dir.create(cur_output_dir)
    }

    cur_config <- config[i,]
    cur_config <- cbind(cur_config, top_num=full_config$top_num)
    cur_platform <- cur_config[["platform_name"]]
    cur_platform_abbr <- cur_config[["platform_abbr"]]
    all_platform_abbr <- c(all_platform_abbr, cur_platform_abbr)
    work_dirs[[cur_platform]] <- work_dir
    print(paste0("platform: ", cur_platform))
    geneset_ids <- gmt2list(cur_config[["gmt_file"]])
    score_list <- scan(cur_config[["score_file"]], what="", sep="\n", quiet=TRUE)
    score_list <- strsplit(score_list, "\t")
    names(score_list) <- sapply(score_list, `[[`, 1)
    geneset_info <- lapply(score_list, `[`, c(-1))
    geneset_info <- lapply(geneset_info, "as.double")

    save(geneset_info, geneset_ids, file=file.path(cur_output_dir, "genesets.RData"))
    top_ret <- topGeneSets(cur_config, geneset_ids, geneset_info, cur_output_dir)
    sc_ret <- topSCGeneSets(cur_config, geneset_ids, geneset_info, cur_output_dir, n_threads)

    # order by name
    geneset_ids <- geneset_ids[order(names(geneset_ids))]
    geneset_info <- geneset_info[order(names(geneset_info))]
    set_platfrom <- rep(cur_platform_abbr, length(geneset_ids))
    set_name <- names(geneset_ids)
    set_size <-  unname(sapply(geneset_ids, length))
    set_score <- unname(sapply(geneset_info, "as.double"))
    all_data <- rbind(all_data, data.frame(platform=set_platfrom, name=set_name, size=set_size, score=set_score, stringsAsFactors=FALSE))
  }
  save(all_data, file=file.path(output_dir, "all_data.RData"))
  # plot
  pdf(NULL)
  theme_set(theme_bw())
  myplot <- ggplot(all_data, aes(x=platform, y=score, color=platform)) +
    scale_color_manual(values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")) +
    geom_jitter(aes(size=size), width = 0.25) + guides(color=FALSE) +
    theme(legend.key = element_rect(color=NA, fill = NA))

  ggsave(file.path(output_dir, "all_data.png"), dpi = 300, myplot)
  ggsave(file.path(output_dir, "all_data.pdf"), myplot)
  dev.off()


  # do AP clustering
  ap_data <- list()
  ap_data$output.dir <- output_dir
  ap_data$md5val <- "sumer"

  genesetIds <- list()
  genesetInfo <- list()
  for(i in seq_len(n_platform)){
    cur_config <- config[i,]
    cur_platform <- cur_config[["platform_name"]]
    cur_platform_abbr <- cur_config[["platform_abbr"]]
    cur_sc_output_dir <- work_dirs[[cur_platform]]
    sc_gs <- read.table(file.path(output_dir, cur_sc_output_dir, "sctopsets.txt"), sep="\t", stringsAsFactors=FALSE, header=TRUE)
    for (row in 1:nrow(sc_gs)) {
      tmpname <- paste(cur_platform_abbr, sc_gs[row, "geneset"], sep="_")
      genesetIds[[tmpname]] <- unlist(strsplit(sc_gs[row, "UserID"], ","))
      genesetInfo[[tmpname]] <- list()
      genesetInfo[[tmpname]][["signP"]] <- sc_gs[row, "score" ]
      genesetInfo[[tmpname]][["leNum"]] <- length(genesetIds[[tmpname]])
      # not used
      genesetInfo[[tmpname]][["PValue"]] <-  NULL
      genesetInfo[[tmpname]][["NES"]] <-  sc_gs[row, "score" ]
      # not used
      genesetInfo[[tmpname]][["FDR"]] <- NULL
      genesetInfo[[tmpname]][["platform"]] <- cur_platform
    }
  }

  ap_data$genesetIds <- genesetIds
  ap_data$genesetInfo <- genesetInfo
  ap_results <- affinityPropagation(ap_data)

  # add results to webpage
  cat("<div style=\"margin-left:10px\" id=\"results\">", file=output_index_file, append=TRUE)
  cat("<h4> Original gene sets</h4>\n", file=output_index_file,append=TRUE,sep="")

  cat("<div id=\"originalsets\">\n", file=output_index_file, append=TRUE)
  original_file_link <- "all_data.pdf"
  original_png_file_link <- "all_data.png"
  image_width <- min(300*n_platform, 1200)
  cat("<a href=\"", original_file_link, "\" download><img id=\"original_data\" width=\"", image_width, "px\" src=\"", original_png_file_link, "\" alt=\"original data\"></a>\n</div>\n", file=output_index_file, append=TRUE, sep="")
  cat("\n</div>\n<h4> Top gene sets (up to ", full_config$top_num, ") by set cover </h4>\n", file=output_index_file,append=TRUE,sep="")
  cat("<div id=\"topsets\" style=\"display: grid; grid-template-columns: 1fr 1fr 1fr;\">\n", file=output_index_file, append=TRUE)
  for(i in seq_len(n_platform)) {
    cur_config <- config[i,]
    cur_platform <- cur_config[["platform_name"]]
    cat("<div id=\"platform",i, "\">\n", file=output_index_file, append=TRUE, sep="")
    cat("<h5 style=\"margin-top:0px;\">", cur_platform, "</h5>\n", file=output_index_file, append=TRUE, sep="")
    cat("<div id=\"sctopsets",i, "\">\n", file=output_index_file, append=TRUE, sep="")
    cur_output_dir <- work_dirs[[cur_platform]]
    sc_file_link <- file.path(cur_output_dir, "sctopsets.pdf")
    sc_png_file_link <- file.path(cur_output_dir, "sctopsets.png")
    sc_txt_link <- file.path(cur_output_dir, "sctopsets.txt")
    cat("<a href=\"", sc_file_link, "\" download><img id=\"", paste0("sctopsetsbarplot", i), "\" width=\"600px\" src=\"", sc_png_file_link, "\" alt=\"set cover top sets barplot\"></a>\n</div>\n", file=output_index_file, append=TRUE, sep="")
    cat("<div id=\"", paste0("sc_info", i), "\" style=\"margin-top:20px;margin-bottom:30px;\"><span><a id=\"",
        paste0("sc_download_data", i), "\" href=\"", sc_txt_link, "\" download class=\"chosen-button\">Download data</a></span></div>\n",
        file=output_index_file, append=TRUE, sep="")
    cat("<a href=\"#top\">Go to top</a>\n", file=output_index_file, append=TRUE, sep="")
    cat("</div>\n", file=output_index_file, append=TRUE, sep="")
  }

  cat("</div><h4>Summary of top gene sets from all platforms </h4>\n", file=output_index_file, append=TRUE)
  cat("<div id='apcy_notes' style='border: 1px solid gray; border-radius: 5px;background-color:#F4F6F6;width:1210px;margin-bottom:30px;padding:5px;'>
      <ul>
      <li>Clustering is performed using <a target='_new' href='https://www.psi.toronto.edu/affinitypropagation/'>affinity propagation</a> algorithm.</li>
      <li>Each node represents a gene set. </li>
      <li>Size of nodes represents the relative number of genes in the set. </li>
      <li>Node shape maps to omics platform (e.g., RNA-seq, Proteomics, Phosphoproteomics). </li>
      <li>Node color maps to value of scores. </li>
      </ul>
</div>\n", file=output_index_file, append=TRUE)
  cat('<div id="cy_toolbar">
      <div id="cy_layout_choice" style="float:left; margin-bottom:0px;">Select a layout style:
      <select class="chzn-layout-select" name="cy_layout_select" style="width:250px;">
      <option value="cose-bilkent">CoSE-Bilkent</option>
      <option value="dagre">Dagre</option>
      </select>
  </div>
  <div id="cy_search" style="margin-left:20px;float:left;">Search:
      <select class="chzn-cy-search-select" data-placeholder="search gene set..." name="cy_search_select" style="width:500px;">
      <option value=""></option>\n', file=output_index_file, append=TRUE)

  ap_result_file <- file.path(output_dir, ap_results[['json']])
  ap_results_df <- fromJSON(ap_result_file)
  all_node_edge_ids <- ap_results_df[["vis"]][["data"]][["id"]]
  pattern <- ""
  for(sd in seq_along(all_platform_abbr)){
    pattern <- paste0(pattern, "^", all_platform_abbr[sd])
    if(sd < length(all_platform_abbr)){
      pattern <- paste0(pattern,"|")
    }
  }
  all_node_ids <- all_node_edge_ids[grepl(pattern, all_node_edge_ids)]
  for(node_id in all_node_ids){
    cat("<option value='", node_id, "'>", gsub("_", " ", node_id),"</option>\n", sep="", file=output_index_file, append=TRUE)
  }

  cat("</select>\n</div>\n</div>\n", file=output_index_file, append=TRUE)
  cat("<div id='apcy_info' style='margin-top:0px;margin-bottom:20px;'>\n", file=output_index_file, append=TRUE)
  cat("<span style='display:inline-block; margin-right:5px;'><a id='acpy_download_data' href='", ap_results[['zip']], "' download class='chosen-button'>Download data</a></span>\n", file=output_index_file, sep="", append=TRUE)
  cat("<select data-placeholder='export network as...' class='chzn-cy-export' style='width:200px;' id='cytscp_exp_sel'>
      <option value=''></option>
      <option value='graphml'>graphml</option>
      <option value='json'>json</option>
      <option value='png'>png</option>
      </select>
      <button class='chosen-button' id='cytscp_exp_btn' style='margin-left:10px;'>Go</button>
      </div>\n", file=output_index_file, sep="", append=TRUE)
  cat("<div id='apcy' data-apjson='", ap_results[['json']],"'></div>\n",  sep="", file=output_index_file, append=TRUE)
  cat("<div class='totop' style='margin-top:10px;margin-bottom:50px;float:left'><a href='#top'>Go to top</a></div></div>", file=output_index_file, append=TRUE)

  # embed json to html
  cat("\n<script>\n", file=output_index_file, append=TRUE)
  cat("var apjson='", file=output_index_file, append=TRUE)

  conn <- file(ap_result_file, open="r")
  # should be one line
  line <- readLines(conn, n=1)
  cat(line, file=output_index_file, append=TRUE)
  cat("'", file=output_index_file, append=TRUE)
  close(conn)
  cat("\n</script>\n", file=output_index_file, append=TRUE)
  cat('<iframe src="https://s3-us-west-2.amazonaws.com/sumer-r-app/foot.html" style="border:0px;height:10px;width:100%"></iframe>', file=output_index_file, append=TRUE)
  cat("</body>\n</html>", file=output_index_file, append=TRUE)
}


#' @importFrom jsonlite fromJSON
check_config <- function(config_file) {
  # maximum of number of platforms we can support
  MAX_N_PLATFORM <- 7
  old_config <- fromJSON(config_file)
  if (! "project" %in% names(old_config)) {
    stop("please provide project description")
  }
  if (! "data" %in% names(old_config)) {
    stop("please provide project data information")
  }
  if (! "top_num" %in% names(old_config)) {
    old_config$top_num <- 50
  }
  config <- old_config$data
  n_platform <- nrow(config)
  if (n_platform > MAX_N_PLATFORM) {
    stop(paste0("currently supports up to ", MAX_N_PLATFORM, " platforms"))
  } else if (n_platform < 1) {
    stop("must provide at least data from 1 platform")
  }
  # check if data file valid
  for (i in seq_len(n_platform)) {
    cur_platform <- config[i, "platform"]
    if(!file.exists(config[i, "gmt_file"])) {
       stop(paste0("gmt_file does not exist for platform ", cur_platform))
    } else if (!file.exists(config[i, "score_file"])){
       stop(paste0("score_file does not exist for platform ", cur_platform))
    }
  }
  # TODO: check platform abbr is unique
  return(old_config)
}
