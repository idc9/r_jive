#' Singluar Value Decomposition.
#'
#' Returns a possibly truncated SVD of a data matrix.
#'
#' Wraps the svd function. Removes the extra singluar values if a truncated svd is computed.
#'
#' @param X Matrix.
#' @param rank Integer. Rank of the desired SVD (optional). If rank==0 returns zeros.
#'
#' @return The SVD of X.
get_svd <- function(X, rank=NULL){
    # SVD <- get_svd(X, rank=2)

    if(is.null(rank)){
        svd(X)
    } else if(rank == 0){
        # TODO: what to do
        decomposition <- list()
        decomposition[['u']] <- matrix(0, ncol=1, nrow=dim(X)[1])
        decomposition[['d']] <- 0
        decomposition[['v']] <- matrix(0, ncol=1, nrow=dim(X)[2])
        decomposition

    } else{
        decomposition <- svd(X, nu=rank, nv=rank)
        decomposition[['d']] <- decomposition[['d']][1:rank]
        decomposition
    }

}



#' Truncates an SVD.
#'
#' Removes columns from the U, D, V matrix computed form an SVD.
#'
#'
#' @param decomposition List. List with entries 'u', 'd', and 'v'from the svd function.
#' @param rank List. List with entries 'u', 'd', and 'v'from the svd function.
#'
#' @return The trucated SVD of X.
truncate_svd <- function(decomposition, rank){

    if(rank==0){
        n <- dim(decomposition[['u']])[1]
        d <- dim(decomposition[['v']])[1]
        decomposition[['u']] <- matrix(0, ncol=1, nrow=n)
        decomposition[['d']] <- 0
        decomposition[['v']] <- matrix(0, ncol=1, nrow=d)
    }else{
        decomposition[['u']] <- decomposition[['u']][, 1:rank, drop=FALSE]
        decomposition[['d']] <- decomposition[['d']][1:rank]
        decomposition[['v']] <- decomposition[['v']][, 1:rank, drop=FALSE]
    }

    decomposition
}


#' Reconstruces the original matrix from its SVD.
#'
#' Computes UDV^T to get the approximate (or full) X matrix.
#'
#' @param decomposition List. List with entries 'u', 'd', and 'v'from the svd function.
#'
#' @return Matrix. The original matrix.
svd_reconstruction <- function(decomposition){

    # decomposition rank -- need to truncated singluar values
    r <- dim(decomposition[['u']])[2]

    decomposition[['u']]  %*%
        diag(decomposition[['d']][1:r], nrow=r, ncol=r) %*%
        t(decomposition[['v']])

}
