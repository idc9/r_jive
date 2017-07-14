context('tests basic ajive output')


test_that("The JIVE decomposition is self consistent",{

    # sample test data
    K <- 2
    blocks <- sample_toy_data(n=200, dx=100, dy=500)

    initial_signal_ranks <- c(2, 2)
    jive_decomp <- ajive(blocks, initial_signal_ranks, full=TRUE)

    # Check SVDs
    for(k in 1:K){
        I <- jive_decomp[[k]][['individual']][['full']]
        U <- jive_decomp[[k]][['individual']][['u']]
        D <- jive_decomp[[k]][['individual']][['d']]
        V <- jive_decomp[[k]][['individual']][['v']]
        I_svd <- U %*% diag(D, ncol=length(D), nrow=length(D)) %*%t (V)
        expect_equal(I, I_svd,
                     info=paste0('SVD of I for k= ', k))

        J <- jive_decomp[[k]][['joint']][['full']]
        U <- jive_decomp[[k]][['joint']][['u']]
        D <- jive_decomp[[k]][['joint']][['d']]
        V <- jive_decomp[[k]][['joint']][['v']]
        J_svd <- U %*% diag(D, ncol=length(D)) %*%t (V)
        expect_equal(J, J_svd,
                     info=paste0('SVD of J for k = ', k))
    }

    # check error matrices
    for(k in 1:K){
        X <- blocks[[k]]
        J <- jive_decomp[[k]][['joint']][['full']]
        I <- jive_decomp[[k]][['individual']][['full']]
        E <- jive_decomp[[k]][['noise']]
        expect_equal(X, J + I + E,
                     info=paste0('X = J + I + E for k= ', k))
    }


    # check projections
    for(k in 1:K){

        U <- jive_decomp[['joint_scores']]
        X <- blocks[[k]]
        J <- jive_decomp[[k]][['joint']][['full']]
        expect_equal(U %*% t(U) %*% X, J,
                     info=paste0('J projection for k= ', k))

        I <- jive_decomp[[k]][['individual']][['full']]
        X_orthog <- (diag(dim(X)[1]) - U %*% t(U)) %*% X
        indiv_rank <- jive_decomp[[k]][['individual']][['rank']]
        I_SVD <- get_svd(X_orthog, indiv_rank)
        expect_equal(I, svd_reconstruction(I_SVD),
                     info=paste0('I projection for k= ', k))
    }

    # TODO: add tests for correct ranks

})
