M1ifun <- function(omega,  Psi){
	sig2b <- omega[3]
      beta1hat <- omega[2]
	sig2deltaif <- beta1hat^2*Psi[1,1] + Psi[2,2] - 2*beta1hat*Psi[1,2]
	#M1 <- (Psi[2,2]*sig2b + beta1hat^2*(Psi[1,1]*Psi[2,2] - Psi[1,2]))/(sig2b + sig2deltaif)
	M1 <- Psi[2,2] - (Psi[2,2] - beta1hat*Psi[1,2])^2/(sig2b + sig2deltaif)
	M1
}


betahatfun <- function(W, Y, Psi){
	Z1 <- mean(W*Y) - Psi[1,2]
	Z2 <- mean(Y)
	Z3 <- mean(W)
	Z4 <- mean(W^2) - Psi[1,1]

	Mzz <- matrix(c(1, Z3, Z3, Z4), nrow = 2, byrow = TRUE)
	betahat <- (solve(Mzz)%*%c(Z2, Z1))[,1]
	betahat0 <- betahat[1]; betahat1 <- betahat[2]

	sig2bhat <- mean((Y - betahat0 - betahat1*W)^2 - Psi[2,2]  - betahat1^2*Psi[1,1] + 2*betahat1*Psi[1,2])*n/(n-2)

	c(betahat, sig2bhat)
}


predfun <- function(omega, Y , W){
	betahat1 <- omega[2]; sig2bhat <- omega[3]
	vi <- Y - omega[1] - omega[2]*W
	gamma2 <- (Psi[2,2] - betahat1*Psi[1,2])/(betahat1^2*Psi[1,1] + Psi[2,2] - 2*betahat1*Psi[1,2]+sig2bhat  )
	thetatilde2 <- Y - gamma2*vi
	thetatilde2
}

#### unequal psi:

M1ifununequal <- function(omega,  Psi11, Psi12, Psi22){
	sig2b <- omega[3]
      beta1hat <- omega[2]
	sig2deltaif <- beta1hat^2*Psi11 + Psi22 - 2*beta1hat*Psi12
	#M1 <- (Psi22*sig2b + beta1hat^2*(Psi11*Psi22 - Psi12))/(sig2b + sig2deltaif)
	M1 <- Psi22 - (Psi22 - beta1hat*Psi12)^2/(sig2b + sig2deltaif)
	M1
}


betahatfununequal <- function(W, Y, Psi11, Psi12, Psi22){
	Z1 <- mean(W*Y) - mean(Psi12)
	Z2 <- mean(Y)
	Z3 <- mean(W)
	Z4 <- mean(W^2) - mean(Psi11)

	Mzz <- matrix(c(1, Z3, Z3, Z4), nrow = 2, byrow = TRUE)
	betahat <- (solve(Mzz)%*%c(Z2, Z1))[,1]
	betahat0 <- betahat[1]; betahat1 <- betahat[2]

	sig2bhat <- mean((Y - betahat0 - betahat1*W)^2 - Psi22  - betahat1^2*Psi11 + 2*betahat1*Psi12)*n/(n-2)

	c(betahat, sig2bhat)
}


predfununequal <- function(omega, Y, W, Psi11, Psi12, Psi22){
	betahat1 <- omega[2]; sig2bhat <- omega[3]
	vi <- Y - omega[1] - omega[2]*W
	gamma2 <- (Psi22 - betahat1*Psi12)/(betahat1^2*Psi11 + Psi22 - 2*betahat1*Psi12 + sig2bhat  )
	thetatilde2 <- Y - gamma2*vi
	thetatilde2
}

ladj <- function(sig2b, Wtemp, Ytemp,  Psi11, Psi12, Psi22, betahat){
	vhat <- Ytemp - betahat[1] - betahat[2]*Wtemp
	sig2deltai <- betahat[2]^2*Psi11 + Psi22 - 2*betahat[2]*Psi12
	vars2b <- sig2b + sig2deltai
	-(sum(-0.5*vhat^2/vars2b) - 0.5*sum(log(vars2b)) + log(sig2b) )
}

ladjML <- function(sig2b, Wtemp, Ytemp,  Psi11, Psi12, Psi22, betahat){
	vhat <- Ytemp - betahat[1] - betahat[2]*Wtemp
	sig2deltai <- betahat[2]^2*Psi11 + Psi22 - 2*betahat[2]*Psi12
	vars2b <- sig2b + sig2deltai
	-(sum(-0.5*vhat^2/vars2b) - 0.5*sum(log(vars2b))   )
}


omegahatfunll <- function(Wtemp, Ytemp, Psi11, Psi12, Psi22){
	Z1 <- mean(Wtemp*Ytemp) - mean(Psi12)
	Z2 <- mean(Ytemp)
	Z3 <- mean(Wtemp)
	Z4 <- mean(Wtemp^2) - mean(Psi11)

	Mzz <- matrix(c(1, Z3, Z3, Z4), nrow = 2, byrow = TRUE)
	betahat <- (solve(Mzz)%*%c(Z2, Z1))[,1]
	betahat0 <- betahat[1]; betahat1 <- betahat[2]

	sig2bhat <- optimize(ladj, c(0, 5), Wtemp, Ytemp, Psi11, Psi12, Psi22, betahat)$minimum

	c(betahat, sig2bhat)

}



omegahatfunML <- function(Wtemp, Ytemp, Psi11, Psi12, Psi22){
	Z1 <- mean(Wtemp*Ytemp) - mean(Psi12)
	Z2 <- mean(Ytemp)
	Z3 <- mean(Wtemp)
	Z4 <- mean(Wtemp^2) - mean(Psi11)

	Mzz <- matrix(c(1, Z3, Z3, Z4), nrow = 2, byrow = TRUE)
	betahat <- (solve(Mzz)%*%c(Z2, Z1))[,1]
	betahat0 <- betahat[1]; betahat1 <- betahat[2]

	sig2bhat <- optimize(ladjML, c(0, 5), Wtemp, Ytemp, Psi11, Psi12, Psi22, betahat)$minimum

	c(betahat, sig2bhat)

}





