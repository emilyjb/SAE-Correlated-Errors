ladjMLMultiCov <- function(sig2b, Wtemp, Ytemp,  designcov , betahat ){
        vhat <- Ytemp - betahat[1] - betahat[2]*Wtemp[,1] - betahat[3]*Wtemp[,2]
        sig2deltai <- t(betahat[-1])%*%designcov[2:3,2:3]%*%betahat[-1] + designcov[1,1] - 2*t(betahat[-1])%*%designcov[1,2:3]
        vars2b <- sig2b + sig2deltai[1,1]
        -(sum(-0.5*vhat^2/vars2b) - 0.5*n*(log(vars2b))   )
}
