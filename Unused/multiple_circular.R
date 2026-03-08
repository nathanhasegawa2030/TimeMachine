#R Circular Test

#Test file by Nathan Hasegawa to understand multiple circular-circular regression.

library(circular)
library(CircStats)
library(ggplot2)
library(ggforce)
library(boot)
library(dplyr)

#This is a shell plotting code for unit circle and donut plots.
#It is written by Nathan Hasegawa.

angle_plot <- function(circData, ccol = "#c9cfd4", pcol = "#007fff") {
  #Generates an angle plot of every angle in circData on the unit circle.
  #circData: an object of type circular that includes the angles to be plotted.
  #Must be in RADIANS.
  #ccol: optional argument for the color of the circle.
  #pcol: optional argument for the color of the points.
  plotx <- cos(circData)
  ploty <- sin(circData)
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), 
                 color = "black", size = 2.5) +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), fill = NA, color = ccol,
                linewidth=4) +
    geom_point(aes(x = plotx, y = ploty), color = pcol, size=3) +
    xlim(-1.1,1.1) +
    ylim(-1.1,1.1) +
    theme_void()
}

donut_plot <- function(estData, obsData, ccol="#c9cfd4", 
                       pcol="#007fff", ncol="#f40000") {
  #Generates a "donut plot" (see Jha & Biswas, 2017) that plots both the
  #residuals and estimated data. In this way, we can view predictive power and
  #accuracy by region.
  #estData: the data estimated by a model.
  #obsData: the actual data.
  #pcol: optional argument for the color used when the residuals are POSITIVE.
  #That is, the true angle is greater than the estimated angle (resid > 0).
  #ncol: optional argument for the color used when the residuals are NEGATIVE.
  
  #Reorder the rows at random. This is done for viewability.
  order <- sample(x=1:length(estData), replace=FALSE)
  estData <- estData[order]
  obsData <- obsData[order]
  
  resid <- sapply(obsData - estData, rescale_angle)
  mag <- 1 + cos(resid)
  plotx <- mag * cos(estData)
  ploty <- mag * sin(estData)
  
  ggplot() +
    geom_circle(aes(x0 = 0, y0 = 0, r = 2), fill = NA, color = ccol,
                linewidth=4) +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), fill = NA, color = ccol,
                linewidth=4) +
    geom_segment(aes(x = 0, y = 0, xend = 2, yend = 0), 
                 color = "black", size = 2.5) +
    geom_point(aes(x = plotx, y = ploty, 
                   color = ifelse(resid > 0, "Positive", "Negative")), 
               size=3) +
    scale_color_manual(values = c("Positive"=pcol,"Negative"=ncol),
                       guide = "none") +
    xlim(-2.1,2.1) +
    ylim(-2.1,2.1) +
    theme_void()
  
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


#This helper function was written by Google Gemini.
rescale_angle <- function(angle) {
  angle <- angle %% (2 * pi) # Wrap to [0, 2*pi)
  angle[angle >= pi] <- angle[angle >= pi] - (2 * pi) # Shift angles > pi to the negative range
  angle
}


#The functions below are written by Nathan Hasegawa.

make_cos_df <- function(x,n) {
  #Makes a data frame of cosine terms for multiples of each predictor, up to n.
  np <- dim(x)[2]
  cos_mult <- function (i) cos(i*x)
  cos_terms <- as.data.frame(lapply(1:n, cos_mult))
  cosnames <- function (i) paste0("cos", as.character((i+np-1) %/% np), "x",
                                  as.character(((i+np-1) %% np) + 1))
  colnames(cos_terms) <- sapply(1:(n*np), cosnames)
  cos_terms
}

make_sin_df <- function(x,n) {
  #Makes a data frame of sine terms for multiples of each predictor, up to n.
  np <- dim(x)[2] #np stands for number of predictors
  sin_mult <- function (i) sin(i*x)
  sin_terms <- as.data.frame(lapply(1:n, sin_mult))
  sinnames <- function (i) paste0("sin", as.character((i+np-1) %/% np), "x",
                                  as.character(((i+np-1) %% np) + 1))
  colnames(sin_terms) <- sapply(1:(n*np), sinnames)
  sin_terms
}


circ_lm <- function(y, x, n) {
  #Performs multiple circular regression using n terms.
  predictors <- data.frame(make_cos_df(x,n), make_sin_df(x,n))
  cosy <- cos(y)
  siny <- sin(y)
  
  #Perform the regression
  cos.lm <- lm(cosy~., cbind(cosy, predictors))
  sin.lm <- lm(siny~., cbind(siny,predictors))
  
  out <- list()
  out$y <- y
  out$x <- x
  out$num_terms <- n
  out$cos.coefficients <- cos.lm$coefficients
  out$sin.coefficients <- sin.lm$coefficients
  out$fitted.values <- atan2(sin.lm$fitted.values, cos.lm$fitted.values)
  out$cos.fitted.values <- cos.lm$fitted.values
  out$sin.fitted.values <- sin.lm$fitted.values
  out$cos.model <- cos.lm
  out$sin.model <- sin.lm
  out$residuals <- y - atan2(sin.lm$fitted.values, cos.lm$fitted.values)
  out$residuals <- sapply(out$residuals, rescale_angle)
  #resid <- (y - out$fitted.values) %% (2*pi)
  #out$residuals <- sapply(resid, function (r) ifelse(r > pi, r - (2*pi) , r))
  #Note that these are circular residuals; i.e. angular error only.
  out
}


#(3) Multiple predictors.

#To start, we consider just 2 uncorrelated predictors. (The correlated
#predictors may pose a problem, as we will see.)
#We will also have no interaction terms at this time.
x1 <- circular(runif(200, 0, 2*pi))
x2 <- circular(runif(200, 0, 2*pi))
y <- atan2(0.15*cos(x1) + 0.25*sin(x2), 0.35*sin(2*x2)) 
+ rvonmises(n=200, mu=circular(0), kappa=1)

#This is not a good test because all the angles are small
#y <- acos(0.1*cos(x1) - 0.2 * cos(2*x1) + 0.05 * sin(x1) + 0.1 * sin(2*x1) + 
#            0.2 * cos(x2) - 0.05*cos(2*x2) - 0.05 * sin(x2) + 0.05 * sin(2*x2))
#     + rvonmises(n=200, mu=circular(0), kappa=10)

#10th order model
x <- data.frame(x1,x2)
  
lm_test <- circ_lm(y,x,2)
angle_plot(lm_test$residuals)
angle_plot(y - lm_test$fitted.values)

#The error actually gets worse as we increase the number of terms past 2. That's
#where I suspect that multicollinearity might become a problem!


#(4) Highly correlated predictors.
#This is a test to see whether regularization may be able to "solve" an
#ill-conditioned problem which includes multicollinearity.

#Generate synthetic data. x1 and x2 are highly correlated; x3 is not.
x1 <- circular(runif(500, 0, 2*pi))
x2 <- x1 + rvonmises(n=500, mu=circular(0), kappa=100)
x3 <- circular(runif(500, 0, 2*pi))
y <- atan2(0.2*cos(x1)-0.5*sin(x3)+0.25*sin(2*x3), 0.8*sin(2*x1)) 
     + rvonmises(n=500, mu=circular(0), kappa=50)
x <- data.frame(x1,x2,x3)

mc_train <- circ_lm(y,x,2)
angle_plot(mc_train$residuals, pcol="#007fff")


xmod <- data.frame(x1,x3)
mc_no_x2 <- circ_lm(y,xmod,2)
angle_plot(mc_no_x2$residuals, pcol="#ffa500")

#See how different the coefficients are
common_ind = c(1,2,4,5,7,8,10,11,13)
mc_train$cos.coefficients[common_ind] - mc_no_x2$cos.coefficients

#The presence of all the nonzero terms, and inaccurate estimation of the
#coefficients, even though y is well-distributed across the circle, shows that
#we are clearly overfitting. Below is a demonstration of that with new data:
x1_test <- circular(runif(500, 0, 2*pi))
x2_test <- x1_test + rvonmises(n=500, mu=circular(0), kappa=100)
x3_test <- circular(runif(500, 0, 2*pi))
y_test <- atan2(0.2*cos(x1_test)-0.5*sin(x3_test)+0.25*sin(2*x3_test), 
                0.8*sin(2*x1_test)) + rvonmises(n=500, mu=circular(0), kappa=50)
x_test <- data.frame(x1_test,x2_test,x3_test)

#Compute the predicted y-values
n = mc_train$num_terms
x_test_full <- data.frame(make_cos_df(x_test, n), make_sin_df(x_test, n))
siny_pred <- predict(mc_train$sin.model, x_test_full)
cosy_pred <- predict(mc_train$cos.model, x_test_full)
y_pred <- atan2(siny_pred, cosy_pred)
angle_plot(y_test - y_pred, pcol="#fd1900")

###TODO: TRY Y_TEST WITH 2 PREDICTORS INSTEAD OF 3 (GET RID OF X_2)

###SEE WHAT HAPPENS IF YOU FIT MULTIVARIATE GAUSSIAN AND THEN TAKE ATAN2 OF ITS
###PREDICTION. THIS MIGHT GET PUSHED TO NEXT WEEK. DO THESE IN PARALLEL (THEY ARE
###SEPARATE TASKS, BUT WE WANT TO COMPARE THEM AT THE END).

###WHAT IF WE JUST TAKE ATAN2 OF THE TIMESIGNATURE RESPONSE?

#I don't want to try regularizing this just yet— want to talk with Braun first
#to see if this is a reasonable approach. But we could try applying the
#regularization directly to the lm portions of this

#(5) A more relevant overfitting test
#Here, we have 50 predictors (200 in the linear models), and only two of them 
#affect the results. We use only first order terms here— higher order harmonics 
#are not expected in the Circadian clock data.

x_mat <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x <- as.data.frame(x_mat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(0.6*cos(x1) - 1.4*cos(x2) + 3.3*sin(x1) + 0.5*sin(x2),
           0.9*cos(x1) - 2.3*cos(x2) + 0.1*sin(x1) - 0.4*sin(x2)) 
           + rvonmises(n=500, mu=circular(0), kappa=50)

#First-order model using all predictors
mc_train <- circ_lm(y,x,1)
angle_plot(mc_train$residuals, pcol="#007fff")


#Try this model on test data
x_mat_test <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x_test <- as.data.frame(x_mat_test)
x1_test <- x_test[,1]
x2_test <- x_test[,2]
y_test <- atan2(0.6*cos(x1_test) - 1.4*cos(x2_test) + 3.3*sin(x1_test) + 0.5*sin(x2_test),
           0.9*cos(x1_test) - 2.3*cos(x2_test) + 0.1*sin(x1_test) - 0.4*sin(x2_test)) + rvonmises(n=500, mu=circular(0), kappa=50)
x_test_full <- data.frame(make_cos_df(x_test, 1), make_sin_df(x_test, 1))
siny_pred <- predict(mc_train$sin.model, x_test_full)
cosy_pred <- predict(mc_train$cos.model, x_test_full)
y_pred <- atan2(siny_pred, cosy_pred)
angle_plot(y_test - y_pred, pcol="#fd1900")
donut_plot(y_pred, y_test)
