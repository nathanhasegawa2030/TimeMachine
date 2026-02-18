#R Circular Test — MGaussian

#Test file by Nathan Hasegawa to implement elastic net regularization into
#multiple circular-circular regression, using the SAME METHODOLOGY AS
#TIMEMACHINE (doing both linear models using the same regularization parameters
#as opposed to doing them separately with possibly disparate regularizations).
library(circular)
library(CircStats)
library(ggplot2)
library(ggforce)
library(boot)
library(glmnet)
library(systemfonts)

source("Desktop/NSH_WI26_Rotation/plotting.R")

#The function below is written by Nathan Hasegawa.

circ.lm.mgaussian <- function(y, x, a=seq(0,1,0.1), Nrep=10, 
                              nl=1000, seed=3330) {
  #Performs multiple circular-circular regression using n terms.
  #y: Vector of circular responses.
  #x: Data frame of circular predictors.
  #n: Number of harmonics to consider (terms up to cos and sin of nx).
  #Default is n=1.
  #a: values of alpha to use in cross-validation for GLMnet.
  #Nrep: number of replicates of 10-fold CV to use for each value of alpha
  #nl: number of lambda values to test for each call to cv.glmnet.
  #seed: optional argument for which seed to set before computing folds. Default
  #is 3330.
  
  #Clean the data
  predictors <- as.matrix(cbind(cos(x), sin(x)))
  cosy <- cos(y)
  siny <- sin(y)
  
  #Initialization
  num_a <- length(a)
  lambda.1se <- as.list(rep(0,num_a))
  cv.al <- numeric(Nrep * nl * num_a)
  cv.la <- numeric(Nrep * nl * num_a)
  cv.rep <- numeric(Nrep * nl * num_a)
  cv.mse <- numeric(Nrep * nl * num_a)
  cv.r2 <- numeric(Nrep * nl * num_a)
  cv.n0 <- numeric(Nrep * nl * num_a)
  index = 1
  
  #Generate folds
  set.seed(seed)
  distn = rep(1:10, ceiling(length(y)/10))
  folds = as.list(rep(0,Nrep))
  for (k in 1:Nrep) {
    folds[[k]] <- sample(distn, length(y), replace=FALSE)
  }
  
  for (i in 1:num_a) {
    #CV loop in alpha. cv.glmnet does it in lambda.
    al <- a[i]
    
    for (j in 1:Nrep) {
      #Fit a model
      model <- cv.glmnet(predictors,as.matrix(cbind(cosy,siny)),keep=T,
                         alpha=al,family="mgaussian",nfolds=10,
                         nlambda=nl,standardize=F,foldid=folds[[j]])
      
      #Find lambda.1se
      ind <- model$index[2]
      lambda.1se[[i]][j] <- model$lambda.1se
      
      #Store the data from this replicate
      nl_used = length(model$lambda)
      cv.al[index:(index+nl_used-1)] = al
      cv.la[index:(index+nl_used-1)] = model$lambda
      cv.rep[index:(index+nl_used-1)] = j
      cv.mse[index:(index+nl_used-1)] = model$cvm
      cv.r2[index:(index+nl_used-1)] = (1 - model$cvm / (var(cosy) + var(siny)))
      cv.n0[index:(index+nl_used-1)] = model$nzero
      index = index + nl_used
    }
    print(paste0("Alpha = ", as.character(al), " complete"))
  }
  
  #Handle the storage of error for all values of alpha and lambda
  z = max(which(cv.al != 0))
  cv.al = cv.al[1:z]
  cv.la = cv.la[1:z]
  cv.rep = cv.rep[1:z]
  cv.mse = cv.mse[1:z]
  cv.r2 = cv.r2[1:z]
  cv.n0 = cv.n0[1:z]
  cv.data <- data.frame(cv.al, cv.la, log10(cv.la), cv.rep, cv.mse, cv.r2, cv.n0)
  colnames(cv.data) <- c("alpha", "lambda", "loglambda", "rep", "mse", "r2", "nonzero")
  
  #Find the optimal alpha and lambda values for the cosine and sine model
  amax <- cv.al[which(cv.r2 == max(cv.r2))[1]]
  aind <- which(a == amax)
  lambdamax <- mean(lambda.1se[[aind]])
  
  #Fit the final models for the cross-validated alpha value.
  #Note that we use all the folds here; no CV.
  model <- glmnet(predictors, as.matrix(cbind(cosy, siny)), keep=T,
                  alpha=amax, standardize=F, family="mgaussian",nlambda = nl)
  
  #Return the output
  out <- list()
  out$y <- y
  out$x <- x
  out$a <- a
  out$Nrep <- Nrep
  out$nl <- nl
  out$amax <- amax
  out$lambdamax <- lambdamax
  out$lambda1se <- lambda.1se
  out$cv.data <- cv.data
  out$cos.coefficients <- coefficients(model, s=lambdamax)$cosy
  out$sin.coefficients <- coefficients(model, s=lambdamax)$siny
  pred <- predict(model, predictors, s=lambdamax)
  pred <- pred[,,1]
  out$fitted.cos <- pred[,1]
  out$fitted.sin <- pred[,2]
  out$fitted.values <- atan2(out$fitted.sin, out$fitted.cos)
  out$residuals <- sapply(y - out$fitted.values, rescale_angle)
  out$model <- model
  out
}

predict.circ.mgaussian <- function(model, newx) {
  #Accepts as input an output a new set of ANGULAR predictors (not cosines and
  #sines) and returns the angular responses predicted by the model. 
  predictors <- as.matrix(cbind(cos(newx), sin(newx)))
  pred <- predict(model$model, predictors, s=model$lambdamax)
  pred <- pred[,,1]
  atan2(pred[,2], pred[,1])
}


