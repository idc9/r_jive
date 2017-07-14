#' Estimate the wedin bound for a data matrix.
#'
#' Esimates the wedin boud for a data matrix with the resampling procedure
#'
#' See page 12 of the AJIVE paper for details.
#'
#' @param X Matrix. The data matrix.
#' @param SVD List. The SVD decomposition of X.
#' @param signal_rank Integer. The estimated signal rank of X.
#'
#' @return The estimated wedin bound.
#'
#' @importFrom stats median
get_wedin_bound <- function(X, SVD, signal_rank){

    # resample for U and V
    U_perp <- SVD[['u']][ , -(1:signal_rank)]
    U_sampled_norms <- wedin_bound_resampling(X=X,
                                              perp_basis=U_perp,
                                              right_vectors=FALSE,
                                              num_samples=1000)

    V_perp <- SVD[['v']][ , -(1:signal_rank)]
    V_sampled_norms <- wedin_bound_resampling(X=X,
                                              perp_basis=V_perp,
                                              right_vectors=TRUE,
                                              num_samples=1000)

    # compute upper bound
    # TODO: which way?
    # EV_estimate <- median(V_sampled_norms)
    # UE_estimate <- median(U_sampled_norms)
    # wedin_bound_est - max(EV_estimate, UE_estimate)/sigma_min
    sigma_min <- SVD[['d']][signal_rank]
    wedin_bound_samples <- mapply(function(u, v)  max(u, v)/sigma_min, U_sampled_norms, V_sampled_norms)
    wedin_bound_est <- median(wedin_bound_samples)

    wedin_bound_est
}


#' Resampling procedure for the wedin bound
#'
#' @param X Matrix. The data matrix.
#' @param perp_basis Matrix. Either U_perp or V_perp: the remaining left/right singluar vectors of X after estimating the signal rank.
#' @param right_vectors Boolean. Right multiplication or left multiplication.
#' @param num_samples Integer. Number of vectors selected for resampling procedure.
wedin_bound_resampling <- function(X, perp_basis, right_vectors, num_samples=1000){

    rank <- dim(perp_basis)[2]
    resampled_norms <- rep(0, num_samples)

    for(s in 1:num_samples){

        sampled_col_index <- sample.int(n=dim(perp_basis)[2],
                                        size=rank,
                                        replace=TRUE)


        perp_resampled <- perp_basis[ , sampled_col_index]

        if(right_vectors){
            resampled_projection <- X %*% perp_resampled
        } else{
            resampled_projection <- t(perp_resampled) %*% X
        }

        # operator L2 norm
        resampled_norms[s] <- norm(resampled_projection,
                                   type='2')
    }

    resampled_norms
}

