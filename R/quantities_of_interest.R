#' Returns one of the full JIVE decomposition matrices.
#'
#' get_block_full returns one of J, I or E.
#'
#' @param ajive_output Output from the ajive function.
#' @param k Integer. Which block.
#' @param type String. type is either 'joint', 'individual' or 'noise'.
#'
#' @return A matrix of dimension equal to the dim(block[[k]]).
#'
#' @export
get_block_full <- function(ajive_output, k, type){

    if(! type  %in% c('joint', 'individual', 'noise')){
        stop('type must be: joint, individual, or noise')
    }

    if(type  == 'noise'){
        ajive_output[[k]][[type]]
    } else{
        ajive_output[[k]][[type]][['full']]
    }

}

#' Returns the common normalized scores.
#'
#' Represents the joint signal among the data blocks.
#'
#' @param ajive_output Output from the ajive function.
#'
#' @return A matrix of dimension n x joint-rank that is a basis of the joint space.
#'
#' @export
get_common_normalized_scores <- function(ajive_output){
    ajive_output[['joint_scores']]
}

#' Returns the AJIVE block scores.
#'
#' Representation of a data block.
#'
#' The scores are a lower dimensional representation of a data block (e.g. latent variables).
#'
#' The scores can either be normalized (i.e. just U) or unnormalized (i.e. UD)/
#'
#' If the SVD of I_k is U, D, V then U are the unformalized individual scores for block k.
#'
#' @param ajive_output Output from the ajive function.
#' @param k Integer. Which block.
#' @param type String. type is either 'joint' or 'individual'.
#' @param normalized Boolean. Whether or not to normalize the scores.
#'
#' @return A matrix of dimension n x rank-(joint/individual-k).
#'
#' @export
get_block_scores <- function(ajive_output, k, type, normalized){
    if(! type  %in% c('joint', 'individual')){
        stop('type must be: joint or individual')
    }

    scores <- ajive_output[[k]][[type]][['u']]

    if(!normalized){
        D <- diag(ajive_output[[k]][[type]][['d']],
                  ncol=dim(scores)[2])

        scores <- scores %*% D
    }

    scores
}

#' AJIVE loadings
#'
#' @param ajive_output Output from the ajive function.
#' @param k Integer. Which block.
#' @param type String. type is either 'joint' or 'individual'.
#'
#' @export
get_block_loadings <- function(ajive_output, k, type){
    if(! type  %in% c('joint', 'individual')){
        stop('type must be: joint or individual')
    }

    ajive_output[[k]][[type]][['v']]
}



#' Joint rank.
#'
#' @param ajive_output Output from the ajive function.
#'
#' @export
get_joint_rank <- function(ajive_output){
    ajive_output[[1]][['joint']][['rank']]
}


#' Individual rank of block k.
#'
#' @param ajive_output Output from the ajive function.
#' @param k Integer. Which block.
#'
#' @export
get_individual_rank <- function(ajive_output, k){
    ajive_output[[k]][['individual']][['rank']]
}

