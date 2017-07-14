
#' Samples toy two block JIVE data.
#'
#' \code{sample_toy_data} samples the two block distribution from AJIVE figure 2.
#'
#' Note AJIVE figure 2 uses n=100, dx=100 and dy=10000.
#'
#' @param n Integer. Number of observations. Must be divisible by 20.
#' @param dx Integer. Dimension of X block. Must be divisible by 2.
#' @param dy Integer. Dimension of Y block. Must be divisible by 10.
#' @param only_observations Boolean. Whether or not to include the true decomposition.
#'
#' @return A list of length 2 with the sampled data.
#'     Each list contains: obs, joint, individual, and noise.
#'
#' @examples
#' blocks <- sample_toy_data()
#' X1 = blocks[[1]]
#' @importFrom stats rnorm
#' @export
sample_toy_data <- function(n=200, dx=100, dy=500, only_observations=TRUE){

    # for exact AJIVE figure 2: n=100, dx= 100, dy=10000

    if(n %% 20 != 0){
        stop('n must be divisible by 20')
    }
    if(dx %% 2 != 0){
        stop('dx must be divisible by 2')
    }
    if(dy %% 10 != 0){
        stop('dy must be divisible by 10')
    }


    X_joint <- cbind(rbind(matrix(1, nrow = n/2, ncol = dx/2),
                           matrix(-1, nrow = n/2, ncol = dx/2)),
                     matrix(0, nrow = n, ncol = dx/2))

    X_indiv <- rbind(matrix(-1, nrow = n/4, ncol = dx),
                     matrix(1, nrow = n/4, ncol = dx),
                     matrix(-1, nrow = n/4, ncol = dx),
                     matrix(1, nrow = n/4, ncol = dx))

    X_noise <- matrix(rnorm(n*dx), nrow = n, ncol = dx)

    X_obs <- 5000 * (X_joint + X_indiv + X_noise)

    Y_joint <- cbind(rbind(matrix(-1, nrow = n/2, ncol = dy/5),
                           matrix(1, nrow = n/2, ncol = dy/5)),
                     matrix(0, nrow = n, ncol = 4*dy/5))

    Y_indiv <- cbind(rbind(matrix(1, nrow = n/5, ncol = dy/2),
                           matrix(-1, nrow = n/5, ncol = dy/2),
                           matrix(-1, nrow = n/5, ncol = dy/2),
                           matrix(1, nrow = n/5, ncol = dy/2),
                           matrix(1, nrow = n/5, ncol = dy/2)),
                     rbind(matrix(1, nrow = n/4, ncol = dy/2),
                           matrix(-1, nrow = n/2, ncol = dy/2),
                           matrix(1, nrow = n/4, ncol = dy/2)))

    Y_noise <- matrix(rnorm(n*dy), nrow = n, ncol = dy)

    Y_obs <- Y_joint + Y_indiv + Y_noise

    if(only_observations){
        output <-  list('1'=X_obs, '2'=Y_obs)
    } else{
        output <- list()
        output[['obs']] <- list('1'=X_obs, '2'=Y_obs)
        output[['decomp']] <- list('1'=list('joint'=list('full'=X_joint),
                                            'individual'=list('full'=X_joint),
                                            'noise'=X_noise),
                                   '2'=list('joint'=list('full'=Y_joint),
                                            'individual'=list('full'=Y_joint),
                                            'noise'=Y_noise))

        # list('1' = list('obs'=X_obs, 'joint' = X_joint, 'individual' = X_indiv, 'noise' = X_noise),
        #     '2' = list('obs'=Y_obs, 'joint' = Y_joint, 'individual' = Y_indiv, 'noise' = Y_noise))

    }
    output
}
