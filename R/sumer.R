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
  if(file.exists(output_index_file)){
    overwrite <- readline(prompt="Output directory not empty, do you want to overwrite existing output? (enter Y or N)")
    if(toupper(overwrite) != 'Y'){
      return(0)
    } else{
      file.remove(output_index_file)
    }
  }
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
  sim <- full_config$similarity
  n_platform <- nrow(config)
  all_platform_abbr <- c()

  # plot the original data
  all_data <- data.frame(platfrom=character(), name=character(), size=integer(), score=numeric())

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
  ap_results <- affinityPropagation(ap_data, sim)

  platform_shape_icon <- c("data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAMAAAAM7l6QAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAALFQTFRFAAAALi4u////KCgocnJyU1NTXFxcX19fYmJiZmZmbGxsZWVldnZ2ZGRkVFRULS0tS0tLX19fcHBwgYGBjY2NlJSUU1NTa2tri4uLsbGx0tLSQEBAY2NjioqKwsLCUFBQbm5up6enT09PcnJyuLi4MzMzMjIyY2Njp6eniYmJampqSkpKi4uLsbGxgICA5+fn8vLy9vb27+/v/v7+////6urq+Pj4+/v76enp7u7u0tLS+nIsWQAAAC90Uk5TAAAAAAAAAAAAAAAAAAAAABBFhbjY5RRlv/D+A0fB+wl07AiE9wICSOvAZRG+77h1H7CsAAAAAWJLR0QCZgt8ZAAAAAlwSFlzAAALEwAACxMBAJqcGAAAAONJREFUKM+Nk9cSgjAQRVcUG6h0EARRsILSIpb//zB1wFGUsJ7HnJlMsnsvQEGLEURJVlRVkSVRYFrwSbuj6cbUjOIkiSNzauhap/22rDWznTQjJVnq2DOLfdnu3F2cSIXTwp13C9tben5Ovsh9b9l72v7K9c/kh7PvrvsAg+Fmm5Ma8u1mOABO211ILZedxgGvO4SCs+chCK80fQ0DEA+EykEE6UbXNwnkiK6jIygxXccKqAldJyqmkcuPzU9DPoaMBRkqv29cSeNCR1gckDBhUcSCDDCuq8H43xI9mFQrOCmP71sujvnFJVioAAAAJXRFWHRkYXRlOmNyZWF0ZQAyMDE5LTAxLTIwVDEwOjI5OjM5LTA2OjAwBp+AHAAAACV0RVh0ZGF0ZTptb2RpZnkAMjAxOS0wMS0yMFQxMDoyOTozOS0wNjowMHfCOKAAAAAASUVORK5CYII=",
                           "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAMAAAAM7l6QAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAAZJQTFRFAAAAZmZme3t7ZGRkZWVlampqaGhoa2trfn5+XFxcFhYWX19fYGBgYmJiXl5eZ2dnbGxsbW1tgICAz8/PaWlpNjY2eHh4LCwsV1dXZmZmXFxckpKSODg4d3d3YWFhqqqqTk5OhISEaWlpvLy8V1dXk5OTDg4OcXFxPj4+VFRUYmJicHBwf39/rq6uaWlpY2NjXl5eaWlpeXl5jIyMoqKiurq6Y2NjcHBwqKioycnJU1NTbW1ttbW1SkpKcXFxv7+/UFBQd3d3ysrKVFRUfn5+WFhYhoaGW1tbjo6OXV1dq6urWFhYqqqqXl5etra2ZWVlwcHBa2trcnJyeXl50dHRoqKie3t7AAAAgoKCxMTElZWVcnJyWVlZOzs7ODg4iIiIt7e3ioqKampqU1NTDQ0NX19fcHBwe3t7Y2NjTExMcHBwampqY2NjT09P19fX/Pz85eXl////8fHxz8/P+fn509PT6Ojo+vr639/f+/v7/f391NTU3d3d5ubm/v7+zMzM1tbW4ODg9fX16enp7e3t1dXVCg2e0QAAAG50Uk5TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAHQn2wSTQ+0OtGP5H9EEhwojSHap7wMWNWCRvuH3Ibvx/QRw9gaD+gyX/RSqHrwrzEXvLOk89FD6ZHqQ/uScAqX7038sAwi59b5kGwENsalLDrACGgcESF/KAAAAAWJLR0Rxrwdc4gAAAAlwSFlzAAALEwAACxMBAJqcGAAAAZtJREFUKM9jYEACjEySkkyMDLgAs5S0tBQzLlkWGdm8PFkZFuyyrGxy8vn58nJs7NiNVlAsKCwsUFTAajwjh5JyIRAoK3Fgcx2zimoRSLpIVQWLdk419eJCMChWV+NEluFi0NDU0tbRLYFIl+jqaGtpagCFgUBP38DQyNjEtLSsvBAKystKTU2MjQwN9PUYzMwtLCuK8gvRQH5RhaWFuRmDlbVNZSFWUGljbcXAbWtnX4VNtsrezpabgYHbwdEJm7STowM3yHE8zi7VmLLVLs48YH/x8rq61aDL1ri58vJC/M3H4u5Riypb6+HOwgcLFx52T686ZNk6L092HqTo8PZBlfbxRo4Wfl8/VMP9fPmRU4J/AKp0gD8rQlZAMLAeVbo+UFAAYbZQUAMkMioqIBHTECSEMJ1ROLgRKNYUEhoWFhrSBGQ2Bgsj3CYSHtFc2BIZFR0TGxsTHRXZUtgcES6CMDwuvjUhMSk5RZSBQTQlOSkxoTU+DmE4f2paekammDiEJy6WmZGelYqQZsvOyZVATl4SuTnZbAwUAwB/NrNNJ6AXQgAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAxOS0wMS0yMFQxMDoyOTozOC0wNjowMKDoi6gAAAAldEVYdGRhdGU6bW9kaWZ5ADIwMTktMDEtMjBUMTA6Mjk6MzgtMDY6MDDRtTMUAAAAAElFTkSuQmCC",
                           "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeCAMAAAAM7l6QAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAAURQTFRFAAAAREREampqUlJSW1tbTExMg4ODVFRUX19fT09PV1dXQEBAaWlp8fHxaGhoZmZmZ2dn////bW1tb29vZWVlWFhYdHR0MTExYGBgioqKzs7OUFBQbW1toqKiWlpafX19vb29REREZGRkkpKSU1NTcnJyq6urXV1dg4ODx8fHS0tLaGhompqaV1dXd3d3tbW1YWFhZWVli4uL0NDQYGBgW1tbjo6OTU1NkpKSdnZ2YmJit7e3Tk5OlZWVeXl5ZGRku7u7UVFRmJiYAAAAe3t7ZmZmv7+/U1NTnJycDg4Ofn5+aGhow8PDVVVVoKCgKioqgYGBampqx8fHXV1dkZGRwMDAwMDAZ2dnY2NjYmJiYmJi5+fn////9/f32NjY/v7+7u7u+vr64eHh8/Pz/Pz89PT02tra9vb23d3d+fn55OTk+/v7r5zijQAAAFt0Uk5TAAAAAAAAAAAAAAAAAAAAAAAAAAAAH6MBP8L+DHDmJqX5BE/SE4PuMrb8CGDdHJX1AkDF/gIp2BbJiUb1FtCRTvga1QGZVfoe2gOhXPsj3wWoY/0o3/37BDxOTeeCqigAAAABYktHRBHitT26AAAACXBIWXMAAAsTAAALEwEAmpwYAAABXklEQVQoz6XS6V8BQRjA8bEq0rJDOuTKkask5L6lXDkjamfDOlL+//dtm4SdfdXv5Xw/My/meQD4S0IcHR8fERKATXpyqjs7052eSHG6ozcYaYRoo0G/I8DdPZP5nEFczLnZtLe7qTKL1fY2QHyDN5vVIlu/KrdfOIZo1dBxYZevHtgHTpd7hNYauV1O7phP4bm8YtFW7NWlR8HhAem99o2RoLHv2ksegBt/IDhBmCbBgP8WhMI0EokOh0AkOhXjaTQClLGZGM9iSqCKJ8Q4EVcBKpl6x+t7KkkBqE7P8TxPqyEAZCaL52yG5P6FyuUZnDL5HMUxvCuwOGYLd/D7z4n7Bxw/3BP8SKhi6UOoH6UixTMsVz6F/Fkp828DUlN9FPJjVUP+DJyq1RfbuqjXqOW2wEZTMDW62YBLJg5b7W1utw6J32XTPnW6zxt1O0/a1abCXv/ldaOXfg+Cf/cFYFfMX1Iw/JsAAAAldEVYdGRhdGU6Y3JlYXRlADIwMTktMDEtMjBUMTA6Mjk6NDAtMDY6MDCZwszWAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDE5LTAxLTIwVDEwOjI5OjQwLTA2OjAw6J90agAAAABJRU5ErkJggg==",
                           "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAWCAMAAADgvdz9AAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAAOpQTFRFAAAAZGRkbW1tZ2dnYmJiY2NjaWlpa2trcXFxfn5+oaGhZWVlSkpKaGhoYWFhnZ2dVFRUhoaGKSkpcXFxx8fHYmJip6enVlZWi4uLPz8/dHR0ZGRkra2tWFhYj4+PR0dHeHh4Z2dns7OzZmZmWlpalJSUTExMfHx8aWlpubm5XFxcmpqaT09PgICAbGxsv7+/X19fn5+fUlJShISEFRUVb29vxcXFY2NjpKSkaGho5ubm////+vr66enpzc3N+/v77e3t09PT/f398PDw2NjY/v7+8/Pz3d3d9vb24uLi9/f3j4+Po6OjoaGh3sU68gAAADp0Uk5TAAAAAAAAAAAAAAAAAXw85hW7AoD9RusaxASMUPAhzQeXWvRaKNUKomX3MNwOrHD6OOMStgJ6/Evnw/hwpNMAAAABYktHRDs5DvRsAAAACXBIWXMAAAsTAAALEwEAmpwYAAAA+UlEQVQoz3XQ5ZLCQBBG0R5GkODusCzubo0Ggu/7v85SQEESJufvrfmqegB0iNvjcROwYKNen89LbRaZ+QPzecDP5JUHQ2HEcCjIZVXYI9EF4iIasQvZdCy+xLtlPCaZ54nkCh9WycTXvHCk0utnXqdTDvM8z2Q3+LLJZkzPeS6/xbdtPmfoxPlTUD9ZLfw6dZ9HSbG0Q51dqUjo56ZyZY8G+0r5fR1xVWuaMWu1qus1L5R644Amh0ZdeV7Hmq0jfjm2mo95BdodlOi07wkE7fZOsnzqdakA1h+cUeo86DMYjsZoYTwawmR60VQp7TKdwOx6+7Nwu87+Aaz1VAtqSkSdAAAAJXRFWHRkYXRlOmNyZWF0ZQAyMDE5LTAxLTIwVDEwOjI5OjQxLTA2OjAwP7XHYgAAACV0RVh0ZGF0ZTptb2RpZnkAMjAxOS0wMS0yMFQxMDoyOTo0MS0wNjowME7of94AAAAASUVORK5CYII=",
                           "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAdCAMAAACKeiw+AAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAAFFQTFRFAAAAWFhYYWFhVFRUZGRkcnJydHR0VlZWdnZ2r6+vcHBwvr6+WlpanZ2dbW1ty8vLgYGBj4+PlZWV19fX3Nzc/Pz8////9/f36Ojo9PT09vb2wCJ8IgAAABN0Uk5TAAAACEZ8hQ2O83b7LNx7/bvc5UDOsoIAAAABYktHRBZ80agZAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAY0lEQVQoz93RSRKAIAwEwGFTUTZFQfj/Q/0AMXf6OpXKBkxLKG2WAaOVAOS62XwPZLutEtiPpww9xw44XwuheocQCykGpJeO34Sz0XE7cXU67hdXzfRmJmf2Zq7G3fz/Y7P6ADbLIAFPjmEKAAAAJXRFWHRkYXRlOmNyZWF0ZQAyMDE5LTAxLTIwVDExOjE1OjEyLTA2OjAwGPWMkQAAACV0RVh0ZGF0ZTptb2RpZnkAMjAxOS0wMS0yMFQxMToxNToxMi0wNjowMGmoNC0AAAAASUVORK5CYII=",
                           "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAeBAMAAADJHrORAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAACpQTFRFAAAAXV1dVVVVbW1tTk5OcnJyubm5cXFxu7u7VFRUubm5+/v7/////Pz8iaLSqAAAAAt0Uk5TAAALjwmD+IT4C/f1byodAAAAAWJLR0QMgbNRYwAAAAlwSFlzAAALEwAACxMBAJqcGAAAAH5JREFUGNNjYAADIWNFBiTA6JoWIoDEF6nY3e6IJO3efWZHiQCS9J4zpxEKQNJnkBSApM8gFECkEQog0nAFMGmYApg0VAFCGqIAIQ1RMBUhDVQQyWB1B4l/djEGH109unno9qG7B8O96P7B8C96eGCEF3p4YoQ3enxgxBciPgEAv7RdD/xsFQAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAxOS0wMS0yMFQxMDoyOTo0MS0wNjowMD+1x2IAAAAldEVYdGRhdGU6bW9kaWZ5ADIwMTktMDEtMjBUMTA6Mjk6NDEtMDY6MDBO6H/eAAAAAElFTkSuQmCC",
                           "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAaCAMAAACXfxyGAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAAQJQTFRFAAAAZWVlZmZmampqZGRkhYWFZ2dncnJy////enp6aWlplJSUdnZ2VVVVcnJyh4eHioqKiYmJcXFxYmJiq6urUFBQhoaGampqv7+/WFhYlZWVIiIidHR0X19fpaWlS0tLgICAZmZmt7e3VVVVjo6Ob29vycnJXFxcnp6eeXl5eXl5XV1doKCgcXFxVlZWkZGRaWlpvLy8T09PhISEYWFhrKysPz8/eXl5W1tbnZ2dcHBwysrKVVVVj4+PaGhourq6VVVVe3t78PDw8fHx5+fn////8vLy0dHR+vr64eHh7u7u9/f32dnZ2tra+Pj4zMzM5ubm/Pz87+/v/f39np6eoaGh72rD8AAAAEJ0Uk5TAAAAAAAAAAAAAAAAAAeQ1tTUkEfvEblo+iPUAoo76QusWfYayXr+Nt+vsDjgfx3OYvkPtkbvBZou33z+G8td9w+7LHUOhgAAAAFiS0dECIbelXoAAAAJcEhZcwAACxMAAAsTAQCanBgAAAEfSURBVCjPdZJnd4JAEEUHXRds0cTee++9N3RtYI35/38lKJEoLPfrPfP27MwDuMPYPxzOTwWn48vOgILB6HLzyxd4t8toUDTyeFfkjZXXg56WMfn8RIXfZ3rGo0BwrdbrYOBv3BAKb4iGTTgkv44j0a1Wb6MR/IiOxXeEwi4ek+Ixm0gSKskEiwGl0nu63qdTCLhMluiQzXCQywt6WsjnoCASXcQCFA/6+lCEUvmoZ4/lEkDlpKdPFWkr1RpPt3ytigCb6w26btTN0lpRs3Wm2XOreb8Za2l3Llp76bQt7KMN3R7l72KvK/fByvQHms0Jgz5jlfuAueFIrUdDDitVHE+u7/Y6GStVBGybzr7fmE1t+L/oaL64/bxwW8zl4V/ItJ8bdLPi4gAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAxOS0wMS0yMFQxMDoyOTozOC0wNjowMKDoi6gAAAAldEVYdGRhdGU6bW9kaWZ5ADIwMTktMDEtMjBUMTA6Mjk6MzgtMDY6MDDRtTMUAAAAAElFTkSuQmCC")
  # add results to webpage
  cat("<div style=\"margin-left:10px\" id=\"results\">", file=output_index_file, append=TRUE)
  cat("<h4> Original gene sets</h4>\n", file=output_index_file,append=TRUE,sep="")

  cat("<div id=\"originalsets\">\n", file=output_index_file, append=TRUE)
  original_file_link <- "all_data.pdf"
  original_png_file_link <- "all_data.png"
  image_width <- min(300*n_platform, 1200)
  cat("<a href=\"", original_file_link, "\" download><img id=\"original_data\" width=\"", image_width, "px\" src=\"", original_png_file_link, "\" alt=\"original data\"></a>\n</div>\n", file=output_index_file, append=TRUE, sep="")
  cat("\n</div>\n<h4> Top gene sets (up to ", full_config$top_num, ") by set cover </h4>\n", file=output_index_file,append=TRUE,sep="")
  cat("<div id=\"topsets\" style=\"overflow:scroll;display: grid; grid-template-columns: 1fr 1fr 1fr;\">\n", file=output_index_file, append=TRUE)
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
  cat("<div id='apcy_notes' style='border: 1px solid gray; border-radius: 5px;background-color:#F4F6F6;width:1610px;margin-bottom:30px;padding:5px;'>
      <ul>
      <li>Clustering is performed using <a target='_new' href='https://www.psi.toronto.edu/affinitypropagation/'>affinity propagation</a> algorithm.</li>
      <li>Each node represents a gene set. </li>
      <li>Size of nodes represents the relative number of genes in the set. </li>
      <li>Node shape maps to omics platform. </li>\n", file=output_index_file, append=TRUE)
  cat('<div class="nodeshape">\n', file=output_index_file, append=TRUE)
  for(i in seq_len(n_platform)) {
    cur_config <- config[i,]
    cur_platform <- cur_config[["platform_name"]]
    cat('<div class="nodeshape-left" style="background: url(', platform_shape_icon[i], ') no-repeat">\n', file=output_index_file, sep='', append=TRUE)
    cat('</div>\n', file=output_index_file, append=TRUE)
    cat('<div class="nodeshape-right">', cur_platform, '</div>\n', file=output_index_file, sep='', append=TRUE)
  }
  cat("</div><li>Node color maps to value of scores. </li>
      </ul>
</div>\n", file=output_index_file, append=TRUE)
  cat('<div class="toolbargrid">
      <div id="cy_layout_choice" style=" margin-bottom:0px;">Select a layout style:
      <select class="chzn-layout-select" name="cy_layout_select" style="width:150px;">
      <option value="cose-bilkent">CoSE-Bilkent</option>
      <option value="dagre">Dagre</option>
      </select>
  </div>
  <div id="cy_search" style="margin-left:20px;">Search:
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

  cat("</select>\n</div>\n", file=output_index_file, append=TRUE)
  cat('<div id="cy_label_choice" style="margin-bottom:0px;">Show/hide label:
<select class="chzn-label-select" name="cy_label_select" style="width:200px;">
<option value="showall">Show all labels</option>
<option value="hide1">Hide first level label</option>
<option value="hide12">Hide first and second level label</option>
</select>
</div>
</div>\n', file=output_index_file, append=TRUE)
  cat("<div id='apcy_info' style='margin:0 0 20px 10px;'>\n", file=output_index_file, append=TRUE)
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
  # currently only support Jaccard or Simpson
  if ("similarity" %in% names(old_config)){
    sim <- old_config$similarity
    if (sim != "Jaccard" && sim != "Simpson" && sim != "Dice"){
      stop("similarity can be: Jaccard or Simpson or Dice")
    }
  } else{
    sim <- "Jaccard"
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
