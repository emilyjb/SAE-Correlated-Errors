

####  MSE of proposed predictor:


M1hat <- M1ifununequal(c(betahat, sig2bhat), Psi11, Psi12, Psi22)

omegahats <- sapply(1:n, function(k){ omegahatfunML(W[-k], Y[-k], Psi11 = Psi11[-k], Psi12 = Psi12[-k], Psi22 = Psi22[-k])})

M1hatrep <- apply(omegahats, 2, M1ifununequal, Psi11 = Psi11, Psi12 = Psi12, Psi22 = Psi22)

M1bias <- apply(M1hatrep,2,mean) - M1hat

estsK <- apply(omegahats  , 2, predfununequal, Y, W, Psi11, Psi12, Psi22)	

M2hat <- apply((estsK - apply(estsK, 1, mean))^2, 1, sum)


M1hatses <- c(M1hatses, M1hat)
M1biases <- c(M1biases, M1bias)
M2hats <-rbind(M2hats, M2hat)
vhatomegahats <- rbind(vhatomegahats, 99*diag(cov(t(omegahats))))

