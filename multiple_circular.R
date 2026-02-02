#R Circular Test

#Test file by Nathan Hasegawa to understand multiple circular-circular regression.

library(circular)
library(CircStats)
library(ggplot2)
library(ggforce)
library(boot)
library(dplyr)

#This is a shell plotting code for unit circle plots.

angle_plot <- function(circData, ccol = "#c9cfd4", pcol = "#007fff") {
  #circData: an object of type circular that includes the angles to be plotted.
  #Must be in RADIANS.
  #ccol: optional argument for the color of the circle.
  #pcol: optional argument for the color of the points.
  plotx <- cos(circData)
  ploty <- sin(circData)
  ggplot() +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), fill = NA, color = ccol,
                linewidth=4) +
    geom_point(aes(x = plotx, y = ploty), color = pcol, size=3) +
    xlim(-1.1,1.1) +
    ylim(-1.1,1.1) +
    theme_classic()
}

#And this helper function (written by Daniel Apley, IEMS) will help with
#cross-validation.

CVInd <- function(n,K) { #n is sample size; K is number of parts; returns K-length list of indices for each part
  m<-floor(n/K) #approximate size of each part
  r<-n-m*K
  I<-sample(n,n) #random reordering of the indices
  Ind<-list() #will be list of indices for all K parts
  length(Ind)<-K
  for (k in 1:K) {
    if (k <= r) kpart <- ((m+1)*(k-1)+1):((m+1)*k)
    else kpart<-((m+1)*r+m*(k-r-1)+1):((m+1)*r+m*(k-r))
    Ind[[k]] <- I[kpart] #indices for kth part of data
  }
  Ind
}


#(3) Multiple predictors.

#To start, we consider just 2 uncorrelated predictors. (The correlated
#predictors may pose a problem, as we will see.)
#We will also have no interaction terms at this time.
x1 <- circular(runif(200, 0, 2*pi))
x2 <- circular(runif(200, 0, 2*pi))
y <- 0.4*cos(2*x1) + 0.7*sin(2*x2) - 0.2*sin(x1) +
  rvonmises(n=200, mu=circular(0), kappa=1)

#First order model
cos_df <- data.frame(cos(y), cos(x1), sin(x1), cos(x2), sin(x2))
sin_df <- data.frame(sin(y), cos(x1), sin(x1), cos(x2), sin(x2))

colnames(cos_df) <- c("cosy", "cosx1", "cosx2", "sinx1", "sinx2")
colnames(sin_df) <- c("siny", "cosx1", "cosx2", "sinx1", "sinx2")
cos.lm1 <- lm(cosy~.,cos_df)
sin.lm1 <- lm(siny~.,sin_df)

yhat1 <- atan2(sin.lm1$fitted.values, cos.lm1$fitted.values)
angle_plot(yhat1 - y)

#Second order model
cos_df <- data.frame(cos(y), cos(x1), cos(2*x1), sin(x1), sin(2*x1),
                     cos(x2), cos(2*x2), sin(x2), sin(2*x2))
sin_df <- data.frame(sin(y), cos(x1), cos(2*x1), sin(x1), sin(2*x1),
                     cos(x2), cos(2*x2), sin(x2), sin(2*x2))

colnames(cos_df) <- c("cosy", "cosx1", "cos2x1", "sinx1", "sin2x1",
                              "cosx2", "cos2x2", "sinx2", "sin2x2")
colnames(sin_df) <- c("siny", "cosx1", "cos2x1", "sinx1", "sin2x1",
                      "cosx2", "cos2x2", "sinx2", "sin2x2")
cos.lm2 <- lm(cosy~.,cos_df)
sin.lm2 <- lm(siny~.,sin_df)

yhat2 <- atan2(sin.lm2$fitted.values, cos.lm2$fitted.values)
angle_plot(yhat2 - y)


#Shell function for nth order model. Note that we can easily define 
#interaction terms by multiplying the other terms.
x <- data.frame(x1,x2)

circ_lm <- function(y, x, n) {
  #Performs multiple circular regression using n terms.
  x <- data.frame(x1,x2)
  np <- dim(x)[2] #np stands for number of predictors
  cos_mult <- function (i) cos(i*x)
  sin_mult <- function (i) sin(i*x)
  cos_terms <- as.data.frame(lapply(1:n, cos_mult))
  sin_terms <- as.data.frame(lapply(1:n, sin_mult))
  
  #Name the columns for interpretability
  cosnames <- function (i) paste0("cos", as.character((i+np-1) %/% np), "x",
                                  as.character(((i+np-1) %% np) + 1))
  sinnames <- function (i) paste0("sin", as.character((i+np-1) %/% np), "x",
                                  as.character(((i+np-1) %% np) + 1))
  colnames(cos_terms) <- sapply(1:(n*np), cosnames)
  colnames(sin_terms) <- sapply(1:(n*np), sinnames)
  
  #Make the regression dataframes
  cos_df <- data.frame(cos(y), cos_terms, sin_terms)
  colnames(cos_df)[1] <- "cosy"
  sin_df <- data.frame(sin(y), cos_terms, sin_terms)
  colnames(sin_df)[1] <- "siny"
  
  #Perform the regression
  cos.lm <- lm(cosy~., cos_df)
  sin.lm <- lm(siny~., sin_df)
  
  out <- list()
  out$y <- y
  out$x <- x
  out$cos.coefficients <- cos.lm$coefficients
  out$sin.coefficients <- sin.lm$coefficients
  out$fitted.values <- atan2(sin.lm$fitted.values, cos.lm$fitted.values)
  out$cos.fitted.values <- cos.lm$fitted.values
  out$sin.fitted.values <- sin.lm$fitted.values
  resid <- (y - out$fitted.values) %% (2*pi)
  out$residuals <- sapply(resid, function (r) ifelse(r > pi, r - (2*pi) , r))
  
  #Note that these are circular residuals; i.e. angular error only.
  
  ###TODO: THE RESIDUALS ARE WRONG. THEY NEED TO BE COMPUTED USING THE MIN
  ###FUNCTION
  
  out
}
  
lm_test <- circ_lm(y,x,49)
angle_plot(lm_test$residuals)
angle_plot(y - lm_test$fitted.values)


#(4) Highly correlated predictors.
#This is a test to see whether regularization may be able to "solve" an
#ill-conditioned problem which includes multicollinearity.

x1 <- circular(runif(200, 0, 2*pi))
x2 <- circular(runif(200, 0, 2*pi))
y <- 0.4*cos(2*x1) + 0.7*sin(2*x2) - 0.2*sin(x1) +
  rvonmises(n=200, mu=circular(0), kappa=1)


