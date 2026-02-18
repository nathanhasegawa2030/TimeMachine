#R Circular Test — Separate Models

#Test file by Nathan Hasegawa to implement elastic net regularization into
#multiple circular-circular regression, using separate regularized linear models 
#to estimate the cosine and sine of the response. 
library(circular)
library(CircStats)
library(ggplot2)
library(ggforce)
library(boot)
library(glmnet)
library(systemfonts)

source("Desktop/NSH_WI26_Rotation/plotting.R")

#The function below is written by Nathan Hasegawa.

circ.lm <- function(y, x, a=seq(0,1,0.1), Nrep=10,
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
  cos.lambda.1se <- as.list(rep(0,num_a))
  sin.lambda.1se <- as.list(rep(0,num_a))
  cos.al <- numeric(Nrep * nl * num_a)
  cos.la <- numeric(Nrep * nl * num_a)
  cos.rep <- numeric(Nrep * nl * num_a)
  cos.mse <- numeric(Nrep * nl * num_a)
  cos.r2 <- numeric(Nrep * nl * num_a)
  cos.n0 <- numeric(Nrep * nl * num_a)
  sin.al <- numeric(Nrep * nl * num_a)
  sin.la <- numeric(Nrep * nl * num_a)
  sin.rep <- numeric(Nrep * nl * num_a)
  sin.mse <- numeric(Nrep * nl * num_a)
  sin.r2 <- numeric(Nrep * nl * num_a)
  sin.n0 <- numeric(Nrep * nl * num_a)
  cindex = 1
  sindex = 1
  
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
      cos.model <- cv.glmnet(predictors, cosy, keep=T,alpha=al,family="gaussian",
                             nfolds=10, foldid=folds[[j]], type.measure="deviance", 
                             nlambda=nl, standardize=F)
      sin.model <- cv.glmnet(predictors, siny, keep=T,alpha=al,family="gaussian",
                             nfolds=10, foldid=folds[[j]], type.measure="deviance", 
                             nlambda=nl, standardize=F)
    
      #Find lambda.1se
      cos.lambda.1se[[i]][j] <- cos.model$lambda.1se
      sin.lambda.1se[[i]][j] <- sin.model$lambda.1se
      
      #Store the data from this replicate
      cos_nl_used = length(cos.model$lambda)
      sin_nl_used = length(sin.model$lambda)
      
      cos.al[cindex:(cindex+cos_nl_used-1)] = al
      cos.la[cindex:(cindex+cos_nl_used-1)] = cos.model$lambda
      cos.rep[cindex:(cindex+cos_nl_used-1)] = j
      cos.mse[cindex:(cindex+cos_nl_used-1)] = cos.model$cvm
      cos.r2[cindex:(cindex+cos_nl_used-1)] = 1 - cos.model$cvm / (var(cosy))
      cos.n0[cindex:(cindex+cos_nl_used-1)] = cos.model$nzero
      cindex = cindex + cos_nl_used
      
      sin.al[sindex:(sindex+sin_nl_used-1)] = al
      sin.la[sindex:(sindex+sin_nl_used-1)] = sin.model$lambda
      sin.rep[sindex:(sindex+sin_nl_used-1)] = j
      sin.mse[sindex:(sindex+sin_nl_used-1)] = sin.model$cvm
      sin.r2[sindex:(sindex+sin_nl_used-1)] = 1 - sin.model$cvm / (var(siny))
      sin.n0[sindex:(sindex+sin_nl_used-1)] = sin.model$nzero
      sindex = sindex + sin_nl_used
    }
    print(paste0("Alpha = ", as.character(al), " complete"))
  }
  #Handle the storage of error for all values of alpha and lambda
  cz = max(which(cos.al != 0))
  sz = max(which(sin.al != 0))
  
  cos.al = cos.al[1:cz]
  cos.la = cos.la[1:cz]
  cos.rep = cos.rep[1:cz]
  cos.mse = cos.mse[1:cz]
  cos.r2 = cos.r2[1:cz]
  cos.n0 = cos.n0[1:cz]
  cos.data <- data.frame(cos.al, cos.la, log10(cos.la), cos.rep, cos.mse, cos.r2, cos.n0)
  colnames(cos.data) <- c("alpha", "lambda", "loglambda", "rep", "mse", "r2", "nonzero")
  
  sin.al = sin.al[1:sz]
  sin.la = sin.la[1:sz]
  sin.rep = sin.rep[1:sz]
  sin.mse = sin.mse[1:sz]
  sin.r2 = sin.r2[1:sz]
  sin.n0 = sin.n0[1:sz]
  sin.data <- data.frame(sin.al, sin.la, log10(sin.la), sin.rep, sin.mse, sin.r2, sin.n0)
  colnames(sin.data) <- c("alpha", "lambda", "loglambda", "rep", "mse", "r2", "nonzero")
  

  #Find the optimal alpha and lambda values for the cosine and sine model
  cos.amax <- cos.al[which(cos.r2 == max(cos.r2))[1]]
  sin.amax <- sin.al[which(sin.r2 == max(sin.r2))[1]]
  cos.aind <- which(a == cos.amax)
  sin.aind <- which(a == sin.amax)
  cos.lambdamax <- mean(cos.lambda.1se[[cos.aind]])
  sin.lambdamax <- mean(sin.lambda.1se[[sin.aind]])
  
  #Fit the final models for the cross-validated alpha value.
  #Note that we use all the folds here; no CV.
  cos.model <- glmnet(predictors,cosy,keep=T,alpha=cos.amax[1], standardize=F,
                      family="gaussian",type.measure="deviance", nlambda=nl)
  sin.model <- glmnet(predictors,siny,keep=T,alpha=sin.amax[1], standardize=F,
                      family="gaussian",type.measure="deviance", nlambda=nl)
  
  #Return the output
  out <- list()
  out$y <- y
  out$x <- x
  out$a <- a
  out$Nrep <- Nrep
  out$nl <- nl
  out$cos.amax <- cos.amax
  out$sin.amax <- sin.amax
  out$cos.lambdamax <- cos.lambdamax
  out$sin.lambdamax <- sin.lambdamax
  out$cos.data <- cos.data
  out$sin.data <- sin.data
  out$cos.coefficients <- coefficients(cos.model, s=cos.lambdamax)
  out$sin.coefficients <- coefficients(sin.model, s=sin.lambdamax)
  out$fitted.cos <- predict(cos.model, predictors, s=cos.lambdamax)
  out$fitted.sin <- predict(sin.model, predictors, s=sin.lambdamax)
  out$fitted.values <- atan2(out$fitted.sin, out$fitted.cos)
  out$residuals <- sapply(y - out$fitted.values, rescale_angle)
  out$cos.model <- cos.model
  out$sin.model <- sin.model
  out
}

predict.circ <- function(model, newx) {
  #Accepts as input an output a new set of ANGULAR predictors (not cosines and
  #sines) and returns the angular responses predicted by the model. 
  predictors <- as.matrix(cbind(cos(newx), sin(newx)))
  predcos <- predict(model$cos.model, predictors, s=model$cos.lambdamax)
  predsin <- predict(model$sin.model, predictors, s=model$sin.lambdamax)
  atan2(predsin, predcos)
}

