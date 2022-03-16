rm(list = ls(all = TRUE))
##### Set the working directory to the appropriate folder:
 
library("sae")
library("saeME")
source("ladjmulticov.R")
source("repestfuns.R")
source("predfunmulticovs.R")

outeriter <- 0

n <- 2000
 
R <- rbind( c(1, 0.2, 0.2), c(0.2, 1, 0.2), c(0.2, 0.2, 1))
Psivec <- sqrt(c(0.25, 0.5, 0.74))
designcov <- diag(Psivec)%*%R%*%diag(Psivec)

beta0 <- 1
beta1 <- 2
beta2 <- 3


sigb <- 1.2/2

x <- rchisq(n, df = 5)
x2 <- rchisq(n, df = 5)

betas <- c()
sig2bhats <- c()
thetatildes <-  c()
thetatilde2s  <-  c()

Ys <- c()
thetas <- c()
estsyl <- c()
thetahatfixeds <- c()

M1hatses<- c()

M1biases <- c()

M2hats <- c()

omegahatsfixed <- c()
mseFHs <- c()
vhatomegahats <- c()
thetahatopts <- c()
predyls <- c()
thetahatkims <- c()

iter <- 0

repeat{


iter <- iter + 1

	svdpsii <- eigen(designcov)
	sqrtpsi <- svdpsii$vectors%*%diag(sqrt(svdpsii$values))%*%t(svdpsii$vectors)
	Y <- rep(0, n)
	W1 <- rep(0, n)
	W2 <- rep(0, n)
	theta <- rep(0, n)
	for(i in 1:n){
		eu <- sqrtpsi%*%rnorm(3)
		b <- rnorm(1, mean = 0, sd = sigb)
		Y[i] <- beta0 + beta1*x[i] + beta2*x2[i] + b + eu[1]
		W1[i] <- x[i] + eu[2]
		W2[i] <- x2[i] + eu[3]
		theta[i] <- beta0 + beta1*x[i] + beta2*x2[i] + b
	}


### Estimate parameters:

Z11 <- mean(W1*Y) - designcov[1,2]
Z12 <- mean(W2*Y) - designcov[1,3]
Z2 <- mean(Y)
Z31 <- mean(W1)
Z32 <- mean(W2)
W <- cbind(W1, W2)
Z4 <- t(W)%*%W/n  - designcov[2:3,2:3]

Mzz <- rbind( c(1, Z31, Z32), cbind(c(Z31, Z32), Z4))
betahat <- (solve(Mzz)%*%c(Z2, Z11, Z12))[,1]

betahat0 <- betahat[1]; betahat1 <- betahat[2];  betahat2 <- betahat[3]

initsig2b <- summary(lm(Y~W))$sigma^2

sig2bhat <- optim(initsig2b, ladjMLMultiCov,    lower = 0, upper = Inf, method = "L-BFGS-B", W=W, Y=Y, designcov = designcov, betahat=betahat)$par

##### sig2bhat <- mean((Y - betahat0 - betahat1*W)^2 - Psi22  - betahat1^2*Psi11 + 2*betahat1*Psi12)*n/(n-2)

sig2bhats <- c(sig2bhats, sig2bhat)
betas <- rbind(betas, betahat )
 

#### Construct predictors:

sig2deltai <- t(betahat[-1])%*%designcov[2:3,2:3]%*%betahat[-1] + designcov[1,1] - 2*t(betahat[-1])%*%designcov[1,2:3]
vi <- Y - betahat0 - betahat1*W[,1] - betahat2*W[,2]

gamma2 <- ((designcov[1,1] - t(betahat[-1])%*%designcov[1,2:3])/(t(betahat[-1])%*%designcov[2:3, 2:3]%*%betahat[-1] + designcov[1,1] - 2*betahat[-1]%*%designcov[1,2:3] +sig2bhat ))[1,1]
thetatilde2 <- Y - gamma2*vi

#### Optimal weighted average of Yi and beta0 + beta1*Wi:

gammaiopt <- ((sig2bhat + t(betahat[-1])%*%designcov[2:3,2:3]%*%betahat[-1] - t(betahat[-1])%*%designcov[1,2:3])/(designcov[1,1] + sig2bhat + t(betahat[-1])%*%designcov[2:3,2:3]%*%betahat[-1] - 2*t(betahat[-1])%*%designcov[1,2:3]))[1,1]
thetahatopt <- gammaiopt*Y + (1 - gammaiopt)*(betahat0 + betahat1*W[,1] + betahat2*W[,2])

#cbind(thetahatopt, thetatilde2)

thetahatopts <- rbind(thetahatopts, thetahatopt)

thetatilde2s <- rbind(thetatilde2s, thetatilde2)
Ys <- rbind(Ys, Y)
thetas <- rbind(thetas, theta)

source("msecorrelatederrorunequalpsiMLMultiCov.R")

print(paste(iter))


if(iter == 1000){break}

}


apply(betas, 2, mean)
mean(sig2bhats)

mean(apply((thetatilde2s - thetas)^2, 2, mean))
mean( M1hatses + M2hats  ) - mean(M1biases)

 

 
