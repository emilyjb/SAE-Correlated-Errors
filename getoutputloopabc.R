
rm(list= ls(all = TRUE))

library("xtable")

outfixpar <- c()
outmse <- c()

#fnames <- c("June2021SimsML/n100UnequalPsi.Rdata","June2021SimsML/n500UnequalPsi.Rdata",
#		"June2021SimsML/n100Cor2.Rdata","June2021SimsML/n500Cor2.Rdata",
#		"June2021SimsML/n100UnequalPsiFlipped.Rdata","June2021SimsML/n500UnequalPsiFlipped.Rdata",
#		"June2021SimsML/n100FlippedCor2.Rdata","June2021SimsML/n500FlippedCor2.Rdata")
#


fnames <- paste("abrhosimsUnequal/Config", 1:8, ".Rdata")


iterfname <- 0

repeat{

iterfname <- iterfname + 1

load(fnames[iterfname])

####  Output for equal psi:

##### Fixed parameter estimates:

tabfixpar <- cbind(c(apply(betas[1:1000,], 2, mean), mean(sig2bhats[1:1000])), c(apply(betas[1:1000,], 2, sd), sd(sig2bhats[1:1000])) ) 
tabfixparyl <- cbind(apply(estsyl[1:1000,],2,mean), apply(estsyl[1:1000,],2,sd))
tabfixparFH <- cbind(apply(omegahatsfixed[1:1000,] , 2, mean), apply(omegahatsfixed[1:1000,] , 2, sd))

outfixpar <- rbind(outfixpar, cbind(n,tabfixpar, tabfixparyl, tabfixparFH) )


#xtable(cbind(tabfixpar, tabfixparyl, tabfixparFH), digits = 3)


###  MC MSE of predictor:

msecompare  <- apply(cbind( apply((Ys[1:1000,] - thetas[1:1000,])^2,2,mean), apply((thetatilde2s[1:1000,] - thetas[1:1000,])^2,2,mean),      apply((predyls - thetas)^2,2,mean), apply( (thetahatfixeds - thetas)^2,2,mean)), 2, mean)

###  Average of estimated MSE: 
meadj <- mean(apply((thetatilde2s - thetas)^2,2,mean))
fh <- mean(apply( (thetahatfixeds - thetas)^2,2,mean))


outmse <- rbind(outmse, c(n,msecompare, mean(M1hatses + M2hats - M1biases) , mean(mseFHs)) )


if(iterfname == length(fnames)){break}

}

a <- c(0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75)
b <- c(0.75, 0.75, 0.75, 0.75, 0.25, 0.25, 0.25, 0.25)
rho <- c(0.2, 0.2, 0.8, 0.8, 0.2, 0.2, 0.8, 0.8)

outfixpar <- cbind(rep(a,each = 3), rep(b, each = 3), rep(rho, each = 3) , outfixpar)
outmse <- cbind(a, b, rho, outmse)


xtable(outfixpar, digits = 3)
xtable(outmse, digits = 3)




