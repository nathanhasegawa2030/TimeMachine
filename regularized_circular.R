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

circ.lm <- function(y, x, a=seq(0,1,0.1), Nrep=3, nl=1000) {
  #Performs multiple circular-circular regression using n terms.
  #y: Vector of circular responses.
  #x: Data frame of circular predictors.
  #n: Number of harmonics to consider (terms up to cos and sin of nx).
  #Default is n=1.
  #a: values of alpha to use in cross-validation for GLMnet.
  #Nrep: number of replicates of 10-fold CV to use for each value of alpha
  #nl: number of lambda values to test for each call to cv.glmnet.
  
  #Clean the data
  predictors <- as.matrix(cbind(cos(x), sin(x)))
  cosy <- cos(y)
  siny <- sin(y)
  
  #Initialization
  num_a <- length(a)
  cos.lambda.1se <- numeric(num_a)
  sin.lambda.1se <- numeric(num_a)
  cos.r2 <- numeric(num_a)
  sin.r2 <- numeric(num_a)
  
  for (i in 1:num_a) {
    #CV loop in alpha. cv.glmnet does it in lambda.
    al <- a[i]
    
    for (j in 1:Nrep) {
      #Fit a model
      cos_model <- cv.glmnet(predictors, cosy, keep=T,alpha=al,family="gaussian",
                             nfolds=10, type.measure="deviance", nlambda=nl,
                             standardize=F)
      sin_model <- cv.glmnet(predictors, siny, keep=T,alpha=al,family="gaussian",
                             nfolds=10, type.measure="deviance", nlambda=nl,
                             standardize=F)
      
      #Find lambda.1se
      cos.ind <- cos_model$index[2]
      sin.ind <- sin_model$index[2]
      cos.lambda.1se[i] <- cos.lambda.1se[i] + (cos_model$lambda[cos.ind])/Nrep;
      sin.lambda.1se[i] <- sin.lambda.1se[i] + (sin_model$lambda[sin.ind])/Nrep;
      
      #Compute CV r^2 at that lambda and update the averages
      cos.r2[i] <- cos.r2[i] + (1 - (cos_model$cvm[cos.ind] / var(cosy)))/Nrep;
      sin.r2[i] <- sin.r2[i] + (1 - (sin_model$cvm[sin.ind] / var(siny)))/Nrep;
    }
    print(paste0("Alpha = ", as.character(al), " complete"))
  }
  
  #Find the optimal alpha values for the cosine and sine model
  cos.amax <- a[which(cos.r2 == max(cos.r2))]
  sin.amax <- a[which(sin.r2 == max(sin.r2))]
  cos.lambdamax <- cos.lambda.1se[which((cos.r2) == max(cos.r2))]
  sin.lambdamax <- sin.lambda.1se[which((sin.r2) == max(sin.r2))]
  
  #Fit the final models for the cross-validated alpha value.
  #Note that we use all the folds here; no CV.
  cos_model <- glmnet(predictors,cosy,keep=T,alpha=cos.amax[1], standardize=F,
                      family="gaussian",type.measure="deviance", nlambda=nl)
  sin_model <- glmnet(predictors,siny,keep=T,alpha=sin.amax[1], standardize=F,
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
  out$cos.r2 <- cos.r2
  out$sin.r2 <- sin.r2
  out$cos.coefficients <- coefficients(cos_model, s=cos.lambdamax)
  out$sin.coefficients <- coefficients(sin_model, s=sin.lambdamax)
  out$fitted.cos <- predict(cos_model, predictors, s=cos.lambdamax)
  out$fitted.sin <- predict(sin_model, predictors, s=sin.lambdamax)
  out$fitted.values <- atan2(out$fitted.sin, out$fitted.cos)
  out$residuals <- sapply(y - out$fitted.values, rescale_angle)
  out$cos.model <- cos_model
  out$sin.model <- sin_model
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

