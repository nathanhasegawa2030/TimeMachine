#R Circular Test

#Test file by Nathan Hasegawa to implement elastic net regularization into
#multiple circular-circular regression.
library(circular)
library(CircStats)
library(ggplot2)
library(ggforce)
library(boot)
library(glmnet)
library(systemfonts)

#This is a shell plotting code for unit circle and donut plots.
#It is written by Nathan Hasegawa.

angle_plot <- function(circData, ccol = "#c9cfd4", pcol = "#007fff", title=NULL) {
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
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust=0.6,
                                    family="Trebuchet MS", size = 20))
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

#This helper function was written by Google Gemini.
rescale_angle <- function(angle) {
  angle <- angle %% (2 * pi) # Wrap to [0, 2*pi)
  angle[angle >= pi] <- angle[angle >= pi] - (2 * pi) # Shift angles > pi to the negative range
  angle
}

#The function below is written by Nathan Hasegawa.

#Note that in this file, unlike in multiple_circular.R, we are assuming that n=1.

circ_lm <- function(y, x, a=seq(0,1,0.1), Nrep=3, nl=1000) {
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


#Generate synthetic data


#(1) Test of what happens when the sine and cosine of the response have
#coefficients with different orders of magnitude, and the response does not
#sample the entire unit circle. This is a great example of
#"garbage in, garbage out."

x_mat <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x <- as.data.frame(x_mat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(6*cos(x1) - 14*cos(x2) + 33*sin(x1) + 4*sin(x2),
           0.2*cos(x1) - 0.21*cos(x2) + 0.01*sin(x1) - 0.4*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)

#We make an angle plot of y. Notice how it is incredibly concentrated towards
#the top and bottom because of the disparate orders of magnitude. As a result,
#this is essentially a classification problem, even though we are trying to
#solve it using regression. The results are predictably bad.
angle_plot(y, pcol="#2e8b56", title="Distribution of Responses")


model <- circ_lm(y,x)

#Plot the residuals
angle_plot(model$residuals, pcol="#007fff", title="Training Residuals")
#donut_plot(yhat, y)

#Generate test data and redo the prediction
x_mat_test <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x_test <- as.data.frame(x_mat_test)
x1_test <- x_test[,1]
x2_test <- x_test[,2]
y_test <- atan2(6*cos(x1_test) - 14*cos(x2_test) + 33*sin(x1_test) + 4*sin(x2_test),
                0.2*cos(x1_test) - 0.21*cos(x2_test) + 0.01*sin(x1_test) - 0.4*sin(x2_test)) + rvonmises(n=500, mu=circular(0), kappa=100)
predictors_test <- as.matrix(cbind(cos(x_test), sin(x_test)))

yhat_test <- atan2(predict(model$sin.model, predictors_test, s=model$sin.lambdamax),
                   predict(model$cos.model, predictors_test, s=model$cos.lambdamax))
resid_test <- sapply(y_test - yhat_test, rescale_angle)
angle_plot(resid_test, pcol="#fd1900", title="Test Residuals")

#For this test data, a better description of the error is a "misclassification
#rate." If the residual is less than pi/2 in absolute value, we say that our
#model correctly classified a test point; otherwise, we say that it didn't.

training_misclass <- sum(abs(model$residuals) > (pi/2)) / length(model$residuals)
test_misclass <- sum(abs(resid_test) > (pi/2)) / length(resid_test)


