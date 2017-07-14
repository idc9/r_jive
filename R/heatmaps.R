#' Heatmap of a single data matrix.
#'
#' @param data Matrix. The data matrix to plot.
#' @param show_color_bar T/F. Whether or not to show the color bar.
#' @param title String. Optional plot title
#' @param xlab String. Optional x axis label.
#' @param ylab String. Optional y axis label.
#'
#'
#' @examples
#' data <- sample_toy_data()
#' Y_obs <- data[[2]][['obs']]
#' data_heatmap(Y_obs)
#'
#' @import ggplot2
#' @export
data_heatmap <- function(data, show_color_bar=TRUE, title='', xlab='', ylab=''){
    # Heatmap of a single data matrix
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The package 'ggplot2' is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop("The package 'reshape2' is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    # could use geom_raster or geom_tile -- I think the former is faster
    # could add to theme:
    # panel.border = element_rect(linetype = 1, size=.5, colour = "grey80")

    reshaped_data <- as.data.frame(reshape2::melt(data))
    colnames(reshaped_data) <- c('obs', 'var', 'value')

    ggplot(data=reshaped_data,
           aes_string(x = 'var', y = 'obs')) +
    geom_raster(aes_string(fill = 'value'), show.legend=show_color_bar) +
    scale_fill_gradient2(low='blue', high='red') +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          legend.position="bottom") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title=title,
         x=xlab,
         y=ylab)

}

#' Heatmaps of several data blocks.
#'
#' @param blocks List. List containing matrices to plot.
#' @param show_color_bar Boolean. Whether or not to display the color bars.
#'
#' @examples
#' library(ggplot2)
#' library(cowplot)
#' data = sample_toy_data()
#' blocks <- lapply(data, function(x) x[['obs']])
#' data_blocks_heatmap(blocks)
#'
#' @export
data_blocks_heatmap <- function(blocks, show_color_bar=TRUE){

    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The package 'cowplot' is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    heatmap_list <- list()
    for(k in 1:length(blocks)){
        heatmap_list[[k]] <- data_heatmap(blocks[[k]],
                                          show_color_bar=show_color_bar,
                                          ylab=ifelse(k==1, 'observation', ''))
    }

    cowplot::plot_grid(plotlist=heatmap_list)
}

#' Heatmaps of JIVE decomposition of several blocks.
#'
#' @param block_decompositions List containing the K JIVE decomposition lists.
#'
#' @examples
#' data <- sample_toy_data()
#' library(cowplot)
#' data_blocks_heatmap(data)
#'
#' @export
decomposition_heatmaps <- function(block_decompositions){
    # plots a heatmap of the given jive decomposotion for each data block

    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The package 'cowplot' is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    K <- length(block_decompositions)

    heatmap_list <- list()
    for(k in 1:K){

        heatmap_list[[k]] <- data_heatmap(block_decompositions[[k]][['obs']],
                                          ylab=ifelse(k==1, 'observations', ''),
                                          show_color_bar=FALSE)

        heatmap_list[[K + k]] <- data_heatmap(block_decompositions[[k]][['joint']],
                                              ylab=ifelse(k==1, 'joint', ''),
                                              show_color_bar=FALSE)

        heatmap_list[[2*K + k]] <- data_heatmap(block_decompositions[[k]][['individual']],
                                                ylab=ifelse(k==1, 'individual', ''),
                                                show_color_bar=FALSE)

        heatmap_list[[3*K + k]] <- data_heatmap(block_decompositions[[k]][['noise']],
                                                ylab=ifelse(k==1, 'noise', ''),
                                                show_color_bar=FALSE)
    }

    cowplot::plot_grid(plotlist=heatmap_list, ncol=K)
}
