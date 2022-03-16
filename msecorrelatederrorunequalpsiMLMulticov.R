

####  MSE of proposed predictor:


M1hat <- M1ifunequalMultiCov(c(betahat, sig2bhat), designcov)

omegahats <- sapply(1:n, function(k){ omegahatfunMLMultiCov(W[-k,], Y[-k], designcov )})

M1hatrep <- apply(omegahats, 2, M1ifunequalMultiCov, designcov)

M1bias <- mean(M1hatrep) - M1hat

estsK <- apply(omegahats  , 2, predfunequalMultiCov, Y, W, designcov)	

M2hat <- apply((estsK - apply(estsK, 1, mean))^2, 1, sum)


M1hatses <- c(M1hatses, M1hat)
M1biases <- c(M1biases, M1bias)
M2hats <-rbind(M2hats, M2hat)
vhatomegahats <- rbind(vhatomegahats, 99*diag(cov(t(omegahats))))

