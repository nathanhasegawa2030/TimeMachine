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

circ.lm.mgaussian <- function(y, x, a=seq(0,1,0.1), Nrep=3, nl=1000) {
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
  lambda.1se <- numeric(num_a)
  r2 <- numeric(num_a)
  
  for (i in 1:num_a) {
    #CV loop in alpha. cv.glmnet does it in lambda.
    al <- a[i]
    
    for (j in 1:Nrep) {
      #Fit a model
      model <- cv.glmnet(predictors,as.matrix(cbind(cosy,siny)),keep=T,
                         alpha=al,family="mgaussian",nfolds=10,
                         nlambda=nl,standardize=F)
      
      #Find lambda.1se
      ind <- model$index[2]
      lambda.1se[i] <- lambda.1se[i] + (model$lambda.1se)/Nrep;
      
      #Compute CV r^2 at that lambda and update the averages
      r2[i] <- r2[i] + (1 - (model$cvm[ind]/(var(cosy) + var(siny))))/Nrep;
    }
    print(paste0("Alpha = ", as.character(al), " complete"))
  }
  
  #Find the optimal alpha and lambda values for the cosine and sine model
  amax <- a[which((r2 == max(r2)))]
  lambdamax <- lambda.1se[which((r2 == max(r2)))]
  
  #Fit the final models for the cross-validated alpha value.
  #Note that we use all the folds here; no CV.
  model <- glmnet(predictors, as.matrix(cbind(cosy, siny)), keep=T,
                  alpha=amax[1], standardize=F, family="mgaussian",nlambda = nl)
  
  #Return the output
  out <- list()
  out$y <- y
  out$x <- x
  out$a <- a
  out$Nrep <- Nrep
  out$nl <- nl
  out$amax <- amax
  out$lambdamax <- lambdamax
  out$r2 <- r2
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


