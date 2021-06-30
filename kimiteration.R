

xmod <- Y
ymod <- W


initlm <- lm(ymod~xmod)
initsig <- summary(initlm)$sigma^2

initmod <- c(coef(initlm), initsig)

curparmod <- initmod

PsiXXmod <- Psi22
PsiYYmod <- Psi11
PsiXYmod <- Psi12

kimiter <- 0

repeat{

	kimiter <- kimiter + 1

	sig2bmod <- curparmod[3]
	beta0mod <- curparmod[1]
	beta1mod <- curparmod[2]

	wbeta1 <- sig2bmod + beta1mod^2*PsiXXmod + PsiYYmod - 2*beta1mod*PsiXYmod

	xbarwmod <- sum( xmod*wbeta1)/sum(wbeta1)
	ybarwmod <- sum( ymod*wbeta1)/sum(wbeta1)

	beta1modt <- sum( wbeta1*((xmod - xbarwmod)*(ymod - ybarwmod) - PsiXYmod))/sum( (xmod - xbarwmod)^2 - PsiXXmod)
	beta0modt <- ybarwmod - xbarwmod*beta1modt

	covtermsige <- beta1mod^2*PsiXXmod + PsiYYmod - 2*beta1mod*PsiXYmod
	kappa1 <- sig2bmod + beta1modt^2*PsiXXmod + PsiYYmod - 2*beta1modt*PsiXYmod

	sig2bmodt <- sum(kappa1*((ymod - beta0mod - beta1mod*xmod)^2  - covtermsige))/sum(kappa1)

	curparmod <- c(beta0modt, beta1modt, sig2bmodt)

	if(kimiter == 1){break}

}

beta1mod <- curparmod[2]
alphakim <- (curparmod[3] + PsiYYmod - beta1mod*PsiXYmod)/( curparmod[3] + PsiYYmod + beta1mod^2*PsiXXmod  - 2*beta1mod*PsiXYmod )
Xtildekim <- (ymod - curparmod[1])/curparmod[2]

thetahatkim <- alphakim*xmod + (1 - alphakim)*Xtildekim

thetahatkims <- rbind(thetahatkims, thetahatkim)








