#' Angle based Joint and Individual Variation Explained
#'
#' Computes the JIVE decomposition.
#'
#' @param blocks List. A list of the data matrices.
#' @param initial_signal_ranks Vector. The initial signal rank estimates.
#' @param full Boolean. Whether or not to store the full J, I, E matrices or just their SVDs (set to FALSE to save memory).
#'
#' @return The JIVE decomposition.
#'
#' @examples
#' data <- sample_toy_data(n=200, dx=100, dy=500)
#' blocks <- lapply(data, function(x) x[['obs']])
#' initial_signal_ranks <- c(2, 2)
#' jive_decomp <- ajive(blocks, initial_signal_ranks)
#'
#' joint_scores <- jive_decomp[['joint_scores']]
#' J_1 <- jive_decomp[[1]][['joint']][['full']]
#' U_individual_2 <- jive_decomp[[2]][['individual']][['u']]
#' individual_rank_2 <- jive_decomp[[2]][['individual']][['rank']]
#'
#' @export
ajive <- function(blocks, initial_signal_ranks, full=TRUE){

    K <- length(blocks)

    # step 1: initial signal space extraction --------------------------------
    # initial estimate of signal space with SVD

    block_svd <- list()
    sv_thresholds <- rep(0, K)
    for(k in 1:K){

        block_svd[[k]] <- get_svd(blocks[[k]])

        sv_thresholds[k] <- get_sv_threshold(singular_values = block_svd[[k]][['d']],
                                               rank=initial_signal_ranks[k])

    }


    # step 2: joint sapce estimation -------------------------------------------------------------

    joint_scores <- get_joint_scores(blocks, block_svd, initial_signal_ranks, sv_thresholds)
    joint_rank <- dim(joint_scores)[2]


    # step 3: final decomposition -----------------------------------------------------

    jive_decomposition <- list()
    for(k in 1:K){
        jive_decomposition[[k]] <- get_final_decomposition(X=blocks[[k]],
                                                           joint_scores=joint_scores,
                                                           sv_threshold=sv_thresholds[k])
    }

    jive_decomposition[['joint_scores']] <- joint_scores
    jive_decomposition[['joint_rank']] <- joint_rank

    jive_decomposition
}

#' The singular value threshold.
#'
#' Computes the singluar value theshold for the data matrix (half way between the rank and rank + 1 singluar value).
#'
#' @param singular_values Numeric. The singular values.
#' @param rank Integer. The rank of the approximation.
get_sv_threshold <- function(singular_values, rank){

    .5 * (singular_values[rank] + singular_values[rank + 1])
}

#' Computes the joint scores.
#'
#' Estimate the joint rank with the wedin bound, compute the signal scores SVD, double check each joint component.
#'
#' @param blocks List. A list of the data matrices.
#' @param block_svd List. The SVD of the data blocks.
#' @param initial_signal_ranks Numeric vector. Initial signal ranks estimates.
#' @param sv_thresholds Numeric vector. The singular value thresholds from the initial signal rank estimates.
#'
#' @return Matrix. The joint scores.
get_joint_scores <- function(blocks, block_svd, initial_signal_ranks, sv_thresholds){

    K <- length(blocks)
    # SVD of the signal scores matrix -----------------------------------------
    signal_scores <- list()
    for(k in 1:K){
        signal_scores[[k]] <- block_svd[[k]][['u']][, 1:initial_signal_ranks[k]]
    }

    M <- do.call(cbind, signal_scores)
    M_svd <- get_svd(M)


    # estimate joint rank with wedin bound -------------------------------------------------------------
    wedin_bounds <- rep(NA, K)
    for(k in 1:K){
        wedin_bounds[k] <- get_wedin_bound(X=blocks[[k]],
                                           SVD=block_svd[[k]],
                                           signal_rank=initial_signal_ranks[k])


    }
    joint_sv_bound <- K - sum(wedin_bounds^2)
    joint_rank_estimate <- sum(M_svd[['d']]^2 > joint_sv_bound)

    # estimate joint score space ------------------------------------

    joint_scores <- M_svd[['u']][ , 1:joint_rank_estimate, drop=FALSE]

    # reconsider joint score space ------------------------------------
    # remove columns of joint_scores that have a
    # trivial projection from one of the data matrices

    to_remove <- c()
    for(k in 1:K){
        for(j in 1:joint_rank_estimate){

            score <- t(blocks[[k]]) %*% joint_scores[ , j]
            sv <- norm(score)

            if(sv < sv_thresholds[[k]]){
                print(paste('removing column', j))
                to_remove <- c(to_remove, j)
                break
            }
        }

    }
    to_keep <- setdiff(1:joint_rank_estimate, to_remove)
    joint_rank <- length(to_keep)
    joint_scores <- joint_scores[ , to_keep, drop=FALSE]

    joint_scores
}

#' Computes the final JIVE decomposition.
#'
#' Computes X = J + I + E for a single data block and the respective SVDs.
#'
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param sv_threshold Numeric vector. The singular value thresholds from the initial signal rank estimates.
#' @param full Boolean. Do we compute the full J, I matrices or just the SVDs (set to FALSE to save memory)..
get_final_decomposition <- function(X, joint_scores, sv_threshold, full=TRUE){

    jive_decomposition <- list()
    jive_decomposition[['individual']] <- get_individual_decomposition(X, joint_scores, sv_threshold, full)
    jive_decomposition[['joint']] <- get_joint_decomposition(X, joint_scores, full)


    if(full){
        jive_decomposition[['noise']] <- X - (jive_decomposition[['joint']][['full']] +
                                                  jive_decomposition[['individual']][['full']])
    } else{
        jive_decomposition[['noise']] <- NA
    }

    jive_decomposition
}

#' Computes the individual matix for a data block.
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param sv_threshold Numeric vector. The singular value thresholds from the initial signal rank estimates.
#' @param full Boolean. Do we compute the full J, I matrices or just the SVD (set to FALSE to save memory).
get_individual_decomposition <- function(X, joint_scores, sv_threshold, full=TRUE){

    X_orthog <- (diag(dim(X)[1]) - joint_scores %*% t(joint_scores)) %*% X

    indiv_decomposition <- get_svd(X_orthog)

    indiv_rank <- sum(indiv_decomposition[['d']] > sv_threshold)

    indiv_decomposition <- truncate_svd(decomposition=indiv_decomposition,
                                        rank=indiv_rank)

    if(full){
        indiv_decomposition[['full']] <- svd_reconstruction(indiv_decomposition)
    } else{
        indiv_decomposition[['full']] <- NA
    }

    indiv_decomposition[['rank']] <- indiv_rank
    indiv_decomposition
}

#' Computes the joint matix for a data block.
#'
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param full Boolean. Do we compute the full J, I matrices or just the SVD (set to FALSE to save memory).
get_joint_decomposition <- function(X, joint_scores, full=TRUE){

    joint_rank <- dim(joint_scores)[2]
    J <-  joint_scores %*% t(joint_scores) %*% X

    joint_decomposition <- get_svd(J, joint_rank)

    if(full){
        joint_decomposition[['full']] <- J
    } else{
        joint_decomposition[['full']] <- NA
    }

    joint_decomposition[['rank']] <- joint_rank
    joint_decomposition

}
