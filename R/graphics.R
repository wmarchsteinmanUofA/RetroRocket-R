plot_kernels <- function(){
  ts <- data.frame(x = 0:8, y = c(0, 1, 1, 0, 1, 1, 0, 1, 0))
  kernel <- data.frame(x = 0:3, y = c(0, 1, 1, 0))
  metrics <- NULL
  plots <- list()
  for (i in 1:((dim(ts)[1] - dim(kernel)[1])+1)){
    k_t <- kernel
    df3 <- ts %>%  mutate(Data = 'Time Series') %>%
      bind_rows(k_t %>%
                  mutate(Data = 'Kernel'))
    df3$Data <- factor(df3$Data, levels = c('Time Series', 'Kernel'))
    plots[[i]] <- ggplot(data = df3, aes(x = x, y = y, col = Data)) + geom_line()
    metrics <- c(metrics, sum(kernel$y * ts[kernel$x + 1,]$y))
    kernel$x <- kernel$x + 1
  }
  ret <- list()
  ret$metrics <- metrics
  ret$plots <- plots
  return(ret)
}
display_kernel_with_ts <- function(kern, ts){
  disp_ts <- pad_ts(ts, kern$padding)

  ggplot() + geom_point()
}
visualization_example <- function(){
  library(pheatmap)
  library(dplyr)
  library(tidyverse)
  library(ComplexHeatmap)
  library(umap)

  foxo_all <- read.csv("~/retrorocket/data/foxo_kernels.csv", header = TRUE)
  p53_all <- read.csv("~/retrorocket/data/p53_kernels.csv", header = TRUE)
  #standardize coefficients
  dilation <- 1
  co_thresh <- 0
  batch <- 2
  run <- 3
  padding <- 0

  #clustering
  l <- 15
  long_p53 <- p53_all[abs(p53_all$coefficients) >= 0.05 & p53_all$dilation == dilation, ]
  long_p53$labels <- NA
  long_p53$labels[long_p53$coefficients < 0] <- "Living"
  long_p53$labels[long_p53$coefficients > 0] <- "Dead"

  long_p53_mat <- long_p53[,c(1:15, 22)]
  for (i in 1:dim(long_p53_mat)[1]){
    r <- unlist(long_p53_mat[i,1:15])
    long_p53_mat[i,1:15] <- (r - mean(r))/sd(r)
  }
  clustered <- kmeans(long_p53_mat, 4, 10)

  c1_names <- names(clustered$cluster[clustered$cluster == 1])
  c2_names <- names(clustered$cluster[clustered$cluster == 2])
  c3_names <- names(clustered$cluster[clustered$cluster == 3])
  c4_names <- names(clustered$cluster[clustered$cluster == 4])
  cluster_hm <- ComplexHeatmap::Heatmap(as.matrix(long_p53_mat[clustered$cluster == 4,1:15]), cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = TRUE, name = "p53_pos_coef", column_title = paste("p53 pos, dila = ",dilation))
  cluster_hm

  cluster_hm <- ComplexHeatmap::Heatmap(as.matrix(long_p53_mat[,2:15]), cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = TRUE, name = "p53_pos_coef", column_title = paste("p53 pos, dila = ",dilation))
  cluster_hm


  u_p53 <- umap(long_p53_mat[,1:15])

  plot.p53(u_p53, factor(clustered$cluster))
  plot.p53(u_p53, factor(long_p53$labels))
  plot.p53(u_p53, factor(long_p53$length))
  plot.p53(u_p53, factor(long_p53$padding))


  plot.p53(u_p53, factor(long_p53[long_p53$labels == "Dead",]$dilation))

  plot.p53 <- function(x, labels,
                       main="Clustered visualizations of p53 kernels",
                       colors=c("#e377c2", "#17becf","#feba32", "#3322bb", "#aa0022","#aabb22","#00bb22","#aa7030"),
                       pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
                       cex.main=1, cex.legend=0.85) {

    layout <- x
    if (is(x, "umap")) {
      layout <- x$layout
    }

    xylim <- range(layout)
    xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
    if (!add) {
      par(mar=c(0.2,0.7,1.2,0.7), ps=10)
      plot(xylim, xylim, type="n", axes=F, frame=F)
      rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
    }
    points(layout[,1], layout[,2], col=colors[as.integer(labels)],
           cex=cex, pch=pch)
    mtext(side=3, main, cex=cex.main)

    labels.u <- unique(labels)
    legend.pos <- "topleft"
    legend.text <- as.character(labels.u)
    if (add) {
      legend.pos <- "bottomleft"
      legend.text <- paste(as.character(labels.u), legend.suffix)
    }

    legend(legend.pos, legend=legend.text, inset=0.03,
           col=colors[as.integer(labels.u)],
           bty="n", pch=pch, cex=cex.legend)
  }

  #end clustering



  p53 <- p53_all[p53_all$dilation == dilation,]

  foxo <- foxo_all[foxo_all$dilation == dilation,]

  colnames(p53)[1:15] <- as.character(0:14*20*dilation)
  colnames(foxo)[1:15] <- as.character(0:14*20*dilation)
  p53_neg <- p53[p53$coefficients < (-co_thresh), ]
  foxo_neg <- foxo[foxo$coefficients < (-co_thresh), ]

  p53_pos <- p53[p53$coefficients > (co_thresh), ]
  foxo_pos <- foxo[p53$coefficients > (co_thresh), ]
  rownames(foxo_neg) <- 1:dim(foxo_neg)[1]
  rownames(p53_neg) <- 1:dim(p53_neg)[1]
  rownames(p53_pos) <- 1:dim(p53_pos)[1]
  rownames(foxo_pos) <- 1:dim(foxo_pos)[1]
  p53_pos_hm <- ComplexHeatmap::Heatmap(as.matrix(p53_pos[,1:15]+p53_pos[,"bias"]), cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE, name = "p53_pos_coef", column_title = paste("p53 pos, dila = ",dilation))
  p53_pos_hm + ComplexHeatmap::Heatmap(as.matrix(foxo_pos[,1:15]+foxo_pos[,"bias"]), cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE, name = "foxo_pos_coef", column_title = paste("foxo pos, dila = ",dilation))



  p53_neg_hm <- ComplexHeatmap::Heatmap(as.matrix(p53_neg[,1:15]+p53_neg[,"bias"]), cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE, name = "p53_neg_coef", column_title = paste("p53 neg, dila = ",dilation))
  p53_neg_hm + ComplexHeatmap::Heatmap(as.matrix(foxo_neg[row_order(p53_neg_hm),1:15]+foxo_neg[,"bias"]), cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, name = "foxo_neg_coef", show_row_names = FALSE, column_title = paste("foxo neg, dila = ",dilation))



  #generate clusters and p53/foxo heatmap plots for all selections based on the following:
  #dilation
  #padding
  #length
  #coeff threshold

  length_list <- unique(p53_all$length)
  dilation_list <- unique(p53_all$dilation)
  padding_list <- unique(p53_all$padding)

  #selection
  for (l in length_list){
    for (d in dilation_list){
      for (p in padding_list){
        selected_p <- p53_all[p53_all$dilation == d & p53_all$length == l & p53_all$padding == p, ]
        selected_f <- foxo_all[foxo_all$dilation == d & foxo_all$length == l & foxo_all$padding == p, ]


      }
    }
  }
  #clustering
  #heatmaps
}

