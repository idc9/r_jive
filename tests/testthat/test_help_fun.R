context('tests helper functions')


test_that("check get_SVD works",{
    X <- matrix(rnorm(n=10*20), nrow=20, ncol=10)
    SVD <- get_svd(X)

    expect_equal(dim(SVD[['u']])[1], dim(X)[1])
    expect_equal(dim(SVD[['u']])[2], min(dim(X)[1], dim(X)[2]))
    expect_equal(dim(SVD[['v']])[1], dim(X)[2])
    expect_equal(dim(SVD[['v']])[2], min(dim(X)[1], dim(X)[2]))

    # make sure lower rank approx works
    # note 1 can give problems
    ranks <- c(1, 5, min(dim(X)))

    for(r in ranks){
        SVD <- get_svd(X, rank=r)
        expect_equal(dim(SVD[['u']])[2], r, info=paste0('r=', r))
        expect_equal(dim(SVD[['u']])[1], dim(X)[1], info=paste0('r=', r))
        expect_equal(length(SVD[['d']]), r, info=paste0('r=', r)) # make sure truncation happend
        expect_equal(dim(SVD[['v']])[2], r, info=paste0('r=', r))
    }


    # check rank 0 works
    SVD <- get_svd(X, rank=0)
    expect_equal(SVD[['u']], matrix(0, ncol=1, nrow=dim(X)[1]))
    expect_equal(SVD[['d']], 0)
    expect_equal(SVD[['v']], matrix(0, ncol=1, nrow=dim(X)[2]))

})

test_that("SVD reconstruction works",{
    X <- matrix(rnorm(n=10*20), nrow=20, ncol=10)
    SVD <- get_svd(X)
    expect_equal(X, svd_reconstruction(SVD))

})


test_that("truncated SVD works",{
    X <- matrix(rnorm(n=10*20), nrow=20, ncol=10)
    SVD <- get_svd(X)

    SVD1 <- truncate_svd(SVD, rank=1) # make sure this doesn't break
    expect_equal(dim(SVD1[['u']])[2], 1)
    expect_equal(length(SVD1[['d']]), 1)
    expect_equal(dim(SVD1[['v']])[2], 1)
    expect_equal(dim(SVD1[['u']])[1], dim(X)[1])
    expect_equal(dim(SVD1[['v']])[1], dim(X)[2])

    SVD5 <- truncate_svd(SVD, rank=5)
    expect_equal(dim(SVD5[['u']])[2], 5)
    expect_equal(length(SVD5[['d']]), 5)
    expect_equal(dim(SVD5[['v']])[2], 5)

    SVD0 <- truncate_svd(SVD, rank=0)
    expect_equal(SVD0[['u']], matrix(0, ncol=1, nrow=dim(X)[1]))
    expect_equal(SVD0[['d']], 0)
    expect_equal(SVD0[['v']], matrix(0, ncol=1, nrow=dim(X)[2]))

})





