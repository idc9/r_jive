context('tests jive output functions')

# sample test data
K <- 2
blocks <- sample_toy_data(n=200, dx=100, dy=500)

initial_signal_ranks <- c(2, 2)
ajive_output <- ajive(blocks, initial_signal_ranks, full=TRUE)

test_that("get_block_full gives proper decomposition",{
    for(k in 1:K){
        X <- blocks[[k]]


        J <- get_block_full(ajive_output, k, type='joint')
        I <- get_block_full(ajive_output, k, type='individual')
        E <- get_block_full(ajive_output, k, type='noise')

        expect_equal(dim(J), dim(X),
                     info=paste0('dim J for k= ', k))
        expect_equal(dim(I), dim(X),
                     info=paste0('dim I for k= ', k))
        expect_equal(dim(E), dim(X),
                     info=paste0('dim E for k= ', k))

        expect_equal(X, J + I + E,
                     info=paste0('X = J + I + E k= ', k))
    }


})


test_that("get_common_normalized_scores",{


    U <- get_common_normalized_scores(ajive_output)

    expect_equal(get_joint_rank(ajive_output), dim(U)[2])

    # check projections
    for(k in 1:K){

        X <- blocks[[k]]
        J <- get_block_full(ajive_output, k, type='joint')

        expect_equal(U %*% t(U) %*% X, J,
                     info=paste0('J projection for k= ', k))


        I <- get_block_full(ajive_output, k, type='individual')
        indiv_rank <- get_individual_rank(ajive_output, k)

        X_orthog <- (diag(dim(X)[1]) - U %*% t(U)) %*% X
        I_SVD <- get_svd(X_orthog, indiv_rank)

        expect_equal(I, svd_reconstruction(I_SVD),
                     info=paste0('I projection for k= ', k))
    }
})


test_that("normalization works for get_block_scores",{
    for(k in 1:K){
        U <- get_block_scores(ajive_output, k, type='joint', normalized=TRUE)
        col_norms <- apply(U, 2, function(c) sqrt(sum(c^2)))

        expect_equal(rep(1, dim(U)[2]),col_norms,
                     info=paste0('U norm k= ', k))

    }
})


test_that("The AJIVE get_ functions are consistent",{


    # Check SVDs
    for(k in 1:K){
        I <- get_block_full(ajive_output, k, type='individual')
        UD <- get_block_scores(ajive_output, k, type='individual', normalized=FALSE)
        V <- get_block_loadings(ajive_output, k, type='individual')

        I_svd <- UD %*%t(V)
        expect_equal(I, I_svd,
                     info=paste0('SVD of I for k= ', k))

        J <- get_block_full(ajive_output, k=k, type='joint')
        UD <- get_block_scores(ajive_output, k, type='joint', normalized=FALSE)
        V <- get_block_loadings(ajive_output, k, type='joint')

        J_svd <- UD %*%t(V)
        expect_equal(J, J_svd,
                     info=paste0('SVD of J for k = ', k))
    }


})
