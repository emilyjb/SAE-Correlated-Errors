rm(list = ls(all = TRUE))
##### Set the working directory to the appropriate folder that contains the repestfuns.R file:
setwd("C:/Users/Emily/Documents/GitHub/SAE-Correlated-Errors")
library("sae")
library("saeME")

source("repestfuns.R")

as <- c(0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75)
bs <- c(0.75, 0.75, 0.75, 0.75, 0.25, 0.25, 0.25, 0.25)
rho <- c(0.2, 0.2, 0.8, 0.8, 0.2, 0.2, 0.8, 0.8)
ns <- rep(c(100, 500), times = 4)

outeriter <- 0

repeat{


outeriter <- outeriter + 1

apick <- as[outeriter]
bpick <- bs[outeriter]
rhopick <- rho[outeriter]
n <- ns[outeriter]


R <- matrix(c(1, rhopick, rhopick, 1), nrow = 2, byrow = TRUE)
D1 <- sqrt(diag(c(apick, bpick)))
D2 <- sqrt(diag(c(apick, bpick)))*0.75
D3 <- sqrt(diag(c(apick, bpick)))*1.25
D4 <- sqrt(diag(c(apick, bpick)))*1.5


Psi1 <- D1%*%R%*%t(D1)
Psi2 <- D2%*%R%*%t(D2)
Psi3 <- D3%*%R%*%t(D3)
Psi4 <- D4%*%R%*%t(D4)

psilist <- vector("list", 4)
psilist[[1]] <- Psi1
psilist[[2]] <- Psi2
psilist[[3]] <- Psi3
psilist[[4]] <- Psi4

	
Psi11 <- rep(c(Psi1[1,1], Psi2[1,1], Psi3[1,1], Psi4[1,1]), each = n/4)
Psi22 <- rep(c(Psi1[2,2], Psi2[2,2], Psi3[2,2], Psi4[2,2]), each = n/4)
Psi12 <- rep(c(Psi1[1,2], Psi2[1,2], Psi3[1,2], Psi4[1,2]), each = n/4)

beta0 <- 1
beta1 <- 2

sigb <- 1.2/2

x <- rchisq(n, df = 5)

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

sqrtpsilist <- vector("list", 4)
uelist <- vector("list", 4)
ue <- c()
for( i in 1:4){ 
	svdpsii <- eigen(psilist[[i]])
	sqrtpsilist[[i]] <- svdpsii$vectors%*%diag(sqrt(svdpsii$values))%*%t(svdpsii$vectors)
      randomgen <- rt(n/2, df = 5)/sqrt(5/3)
	randomgen <- matrix(randomgen, nrow = 2)
	uelist[[i]] <- sqrtpsilist[[i]]%*%randomgen
	ue <- cbind(ue, uelist[[i]])
}
u <- ue[1,]
e <- ue[2,]

b <- rt(n, df = 5)/sqrt(5/3)*sigb  

Y <- beta0 + beta1*x + b + e
W <- x + u

### Estimate parameters:

Z1 <- mean(W*Y) - mean( Psi12)
Z2 <- mean(Y)
Z3 <- mean(W)
Z4 <- mean(W^2) - mean(Psi11)

Mzz <- matrix(c(1, Z3, Z3, Z4), nrow = 2, byrow = TRUE)
betahat <- (solve(Mzz)%*%c(Z2, Z1))[,1]

betahat0 <- betahat[1]; betahat1 <- betahat[2]

sig2bhat <- optimize(ladjML, c(0, 5), W, Y, Psi11, Psi12, Psi22, betahat)$minimum


##### sig2bhat <- mean((Y - betahat0 - betahat1*W)^2 - Psi22  - betahat1^2*Psi11 + 2*betahat1*Psi12)*n/(n-2)

sig2bhats <- c(sig2bhats, sig2bhat)
betas <- rbind(betas, betahat )
 

#### Construct predictors:

sig2deltai <- betahat1^2*Psi11 + Psi22 - 2*betahat1*Psi12
vi <- Y - betahat0 - betahat1*W

gamma2 <- (Psi22 - betahat1*Psi12)/(betahat1^2*Psi11 + Psi22 - 2*betahat1*Psi12 +sig2bhat )
thetatilde2 <- Y - gamma2*vi


#### Optimal weighted average of Yi and beta0 + beta1*Wi:

gammaiopt <- (sig2bhat + betahat[2]^2*Psi11 - betahat[2]*Psi12)/(Psi22 + sig2bhat + betahat[2]^2*Psi11 - 2*betahat[2]*Psi12)
thetahatopt <- gammaiopt*Y + (1 - gammaiopt)*(betahat0 + betahat1*W)


thetahatopts <- rbind(thetahatopts, thetahatopt)

thetatilde2s <- rbind(thetatilde2s, thetatilde2)
Ys <- rbind(Ys, Y)
thetas <- rbind(thetas, beta0 + beta1*x + b)

print(paste(iter))


#####  Ybarrah and Lohr Estimator 

fityl <- FHme(Y~W, vardir = Psi22, var.x = Psi11)
predyl <- fityl$eblup

betahatyl <- fityl$fit$estcoef[,"beta"]

predyls <- rbind(predyls, predyl)
estsyl <- rbind(estsyl, c(betahatyl, fityl$fit$refvar))

##### FH treating covariate as fixed:
dfFH <- data.frame(Y, W, vardir = Psi22)
FHfixed <- eblupFH(dfFH$Y~dfFH$W, vardir = dfFH$vardir)
thetahatfixed <- FHfixed$eblup[,1]

betahatsfixed <- FHfixed$fit$estcoef[,1]
sig2bhatfixed <- FHfixed$fit$refvar

thetahatfixeds <- rbind(thetahatfixeds, thetahatfixed)
omegahatsfixed <- rbind(omegahatsfixed, c(betahatsfixed, sig2bhatfixed))
mseFHs <- cbind(mseFHs, mseFH(dfFH$Y~dfFH$W, vardir = dfFH$vardir)$mse)

source("msecorrelatederrorunequalpsiML.R")
 
if(iter == 1000){break}

}

print(paste("config", outeriter))

save.image(paste("abrhosimsUnequalT/Config", outeriter,".Rdata") )

if(outeriter == 8){break}

}

