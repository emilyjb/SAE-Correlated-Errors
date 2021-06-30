
rm(list= ls(all = TRUE))

library("xtable")

outfixpar <- c()
outmse <- c()


fnames <- c("June2021SimsML/n100EqualPsiFlipped.Rdata", "June2021SimsML/n100UnequalPsiFlipped.Rdata","June2021SimsML/n500EqualPsiFlipped.Rdata",  "June2021SimsML/n500UnequalPsiFlipped.Rdata")

iterfname <- 0

repeat{

iterfname <- iterfname + 1

load(fnames[iterfname])

####  Output for equal psi:

##### Fixed parameter estimates:

tabfixpar <- cbind(c(apply(betas, 2, mean), mean(sig2bhats)), c(apply(betas, 2, sd), sd(sig2bhats)) ) 
tabfixparyl <- cbind(apply(estsyl,2,mean), apply(estsyl,2,sd))
tabfixparFH <- cbind(apply(omegahatsfixed , 2, mean), apply(omegahatsfixed , 2, sd))

outfixpar <- rbind(outfixpar, cbind(tabfixpar, tabfixparyl, tabfixparFH) )


#xtable(cbind(tabfixpar, tabfixparyl, tabfixparFH), digits = 3)


###  MC MSE of predictor:

msecompare  <- apply(cbind( apply((Ys - thetas)^2,2,mean), apply((thetatilde2s - thetas)^2,2,mean),  apply((thetahatkims - thetas)^2,2,mean),     apply((predyl - thetas)^2,2,mean), apply( (thetahatfixeds - thetas)^2,2,mean)), 2, mean)

###  Average of estimated MSE: 
meadj <- mean(apply((thetatilde2s - thetas)^2,2,mean))
fh <- mean(apply( (thetahatfixeds - thetas)^2,2,mean))


outmse <- rbind(outmse, c(msecompare, mean(M1hatses + M2hats - M1biases) , mean(mseFHs)) )


if(iterfname == length(fnames)){break}

}






