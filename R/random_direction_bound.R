#' Estimate the wedin bound for a data matrix.
#'
#' Samples from the random direction bound. Returns on the scale of squared singular value.
#'
#' @param n_obs. The number of observations.
#' @param n_features. The number of features in each data matrix
#' @param num_samples Integer. Number of vectors selected for resampling procedure.
#'
#' @return rand_dir_samples
#'
get_random_direction_bound <- function(n_obs, dims, num_samples=1000){
    n_blocks <- length(dims)
    rand_dir_samples <- rep(0, num_samples)
    for(s in 1:num_samples){
        rand_subspaces <- list()
        for(b in 1:n_blocks){
            X <- matrix( rnorm(n_obs * dims[b], mean=0,sd=1), n_obs, dims[b])
            U <- get_svd(X)[['u']]

            rand_subspaces[[b]] <- U

        }
        M <- do.call(cbind, rand_subspaces)
        M_svd <- get_svd(M, rank=min(dims))

        rand_dir_samples[s] <- M_svd[['d']][1]^2

    }

    rand_dir_samples
}



