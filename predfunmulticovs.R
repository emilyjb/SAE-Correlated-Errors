 M1ifunequalMultiCov <- function(omega,  designcov){
      sig2b <- omega[length(omega)]
      beta1hat <- omega[2]
	beta2hat <- omega[3]
	betahat <- omega[-length(omega)]
      sig2deltaif <- t(betahat[-1])%*%designcov[2:3,2:3]%*%betahat[-1] + designcov[1,1] - 2*t(betahat[-1])%*%designcov[1,2:3]
      #M1 <- (Psi22*sig2b + beta1hat^2*(Psi11*Psi22 - Psi12))/(sig2b + sig2deltaif)
      M1 <- designcov[1,1] - (designcov[1,1] - t(betahat[-1])%*%designcov[1,2:3])^2/(sig2b + sig2deltaif[1,1])
      M1
}

omegahatfunMLMultiCov <-  function(Wtemp, Ytemp, designcov ){
	W1 <- Wtemp[,1]; W2 <- Wtemp[,2]; Y <- Ytemp
	Z11 <- mean(W1*Y) - designcov[1,2]
	Z12 <- mean(W2*Y) - designcov[1,3]
	Z2 <- mean(Y)
	Z31 <- mean(W1)
	Z32 <- mean(W2)
	W <- cbind(W1, W2)
	Z4 <- t(W)%*%W/n  - designcov[2:3,2:3]

	Mzz <- rbind( c(1, Z31, Z32), cbind(c(Z31, Z32), Z4))
	betahat <- (solve(Mzz)%*%c(Z2, Z11, Z12))[,1]

      betahat0 <- betahat[1]; betahat1 <- betahat[2]

      sig2bhat <- optimize(ladjMLMultiCov, c(0, 5), Wtemp, Ytemp, designcov, betahat)$minimum

        c(betahat, sig2bhat)

}

predfunequalMultiCov <- function(omega, Y, W, designcov ){

	betahat <- omega[-length(omega)]
	sig2bhat <- omega[length(omega)]
	sig2deltai <- t(betahat[-1])%*%designcov[2:3,2:3]%*%betahat[-1] + designcov[1,1] - 2*t(betahat[-1])%*%designcov[1,2:3]
	vi <- Y - betahat0 - betahat1*W[,1] - betahat2*W[,2]
	gamma2 <- ((designcov[1,1] - t(betahat[-1])%*%designcov[1,2:3])/(t(betahat[-1])%*%designcov[2:3, 2:3]%*%betahat[-1] + designcov[1,1] - 2*betahat[-1]%*%designcov[1,2:3] +sig2bhat ))[1,1]
	thetatilde2 <- Y - gamma2*vi
	thetatilde2 

}




