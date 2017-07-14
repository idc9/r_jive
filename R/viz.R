#' Scree plot
#'
#' @param singular_values Numeric. The singluar values.
#' @param title String. The title of the plot.
#'
#' @examples
#' X <- matrix(rnorm(n=10*20), nrow=20, ncol=10) + matrix(1, nrow=20, ncol=10)
#' singluar_values <- svd(X)[['d']]
#' scree_plot(singluar_values)
#'
#' @export
scree_plot <- function(singular_values, title=''){

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The package 'ggplot2' is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    index <- 1:length(singular_values)
    ggplot(as.data.frame(x=cbind(index, singular_values)),
           aes(x=index, y=singular_values)) +
        geom_point() +
        geom_line() +
        labs(y='singluar value', x='index', title=title)
}


#' Shows the scree plot for each data block.
#'
#' @param blocks The data blocks.
#'
#' @export
scree_plot_blocks <- function(blocks){
    plot_blocks_grid(blocks, plot_fun = function(X, title) scree_plot(get_svd(X)[['d']], title))
}


#' Makes a plot for each block in a hotizontal grid.
#'
#' @param blocks List.
#' @param plot_fun Function. A plotting function that takes a sigle data matrix and makes a plot
#' @param ... Additional arguments to be passed to plot_fun
plot_blocks_grid <- function(blocks, plot_fun, ...){

    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The package 'cowplot' is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    plotslist <- list()
    for(k in 1:length(blocks)){
        plotslist[[k]] <- plot_fun(blocks[[k]], title=paste0('block ', k), ...)
    }
    cowplot::plot_grid(plotlist=plotslist)
}
