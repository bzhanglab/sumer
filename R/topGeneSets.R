#' topGeneSets
#' @import ggplot2
topGeneSets <- function(cur_config, geneset_ids, geneset_info, output_dir, prefix=""){

  # score is signed logP

	pdf(NULL)
	plot.df <- data.frame(name=character(), score=double(), myclass=character(),
												stringsAsFactors=FALSE)
  for(i in seq_along(geneset_info)){
    pos.neg <- ifelse(geneset_info[[i]] >= 0, "a", "b")
    plot.df <- rbind(plot.df, list(name=as.character(names(geneset_info[i])),
                                   score=geneset_info[[i]], myclass=pos.neg), stringsAsFactors=FALSE)
  }

	# select top 50
	max_size <- cur_config[["top_num"]]
	num.rows <- nrow(plot.df)
	if(nrow(plot.df) > max_size){
		plot.df <- plot.df[order(-abs(plot.df["score"])),]
		plot.df <- plot.df[1:max_size,]
		num.rows <- max_size
	}


  topsets.out.data <- data.frame(geneset=character(), score=double(), UserID=character(),stringsAsFactors=FALSE)

	for(name in names(geneset_info)){
		if(name %in% plot.df[,c("name")]){
			cur.val <- geneset_info[[name]]
			tmp.df <- data.frame(geneset=name,
													 score=cur.val,  UserID=paste(geneset_ids[[name]],collapse=","),
													 stringsAsFactors=FALSE)
			topsets.out.data <- rbind(topsets.out.data, tmp.df)
		}
	}
  write.table(topsets.out.data,
              file=file.path(output_dir, paste0(prefix, "topsets.txt")),
              row.names=F,col.names=T,sep="\t",quote=F)

	max.y <- ceiling(max(plot.df[["score"]]))+1
	min.y <- floor(min(plot.df[["score"]]))-1

  if(min.y > 0){
     min.y <- -max.y
  } else if(max.y < 0){
     max.y <- -min.y
  } else{
    maxabs <- max(abs(min.y), abs(max.y))
    min.y <- -maxabs
    max.y <- maxabs
  }

  # geneset name has category information
  # remove that before plotting
  short.name <- unname(sapply(plot.df$name,function(x){strsplit(x,"__",fixed=TRUE)[[1]][2]}))
  plot.df$short.name <- short.name

  `%>%` <- magrittr::`%>%`
	plot.df <- plot.df %>%
		dplyr::mutate(
					 name = factor(name, levels = name[order(score, decreasing = FALSE)]),
					 label_y = ifelse(score < 0, 0.2, -0.2),
					 label_hjust = ifelse(score < 0, 0, 1)
					 )

	my_plot <- ggplot(plot.df, aes(x = name, y = score, fill = myclass)) +
		geom_bar(stat = "identity", col = "black") +
		geom_text(aes(y = label_y, label = short.name, hjust = label_hjust)) +
		coord_flip() +
		scale_fill_manual(values = c(a = "#e78ac3", b = "#8da0cb")) +
		#scale_fill_manual(values = c(a = "forestgreen", b = "goldenrod")) +
		theme_minimal() +
		theme(axis.text.y = element_blank(),
					axis.ticks.y = element_blank(),
					axis.title.y = element_blank(),
					legend.position = "none",
					legend.justification = 0.05,
					legend.title = element_blank(),
					panel.grid.major.y = element_blank(),
					panel.grid.minor.y = element_blank(),
					panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
					panel.grid.minor.x = element_blank()) +
 				 scale_y_continuous(expression(italic("score")),
														breaks = min.y:max.y, limits = c(min.y, max.y))

				 print(my_plot)
				 topset.png <- file.path(output_dir, paste0(prefix,"topsets.png"))
				 topset.pdf <- file.path(output_dir, paste0(prefix,"topsets.pdf"))

				 ggsave(topset.png, width = 10, height = 0.2*num.rows + 1, dpi = 300, my_plot)
				 ggsave(topset.pdf, width = 10, height = 0.2*num.rows + 1, dpi = 300, my_plot)

         dev.off()
				 return(list(topset.num=num.rows, topset.png=basename(topset.png),
										 topset.pdf=basename(topset.pdf)))
}


#' @importFrom uuid UUIDgenerate
uuid.gen <- function(){
  uuid <- UUIDgenerate()
  return(gsub("-","",uuid))
}
