#R Circular Test

#Test file by Nathan Hasegawa to understand how circular regression works.

library(circular)
library(CircStats)
library(ggplot2)
library(ggforce)
library(boot)

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

#(1) A simple case.
#In this first example, cos(y) = 0.4cos(2x) + 0.7sin(2x).
#The sine term is irrelevant.
x <- circular(runif(200, 0, 2*pi))
y <- acos(0.4 * cos(2*x) + 0.7 * sin(2*x)) +
  rvonmises(n=50, mu=circular(0), kappa=100)

# Fit a circular-circular regression model.
circ.lm <- lm.circular(y, x, type="c-c", order=2)
yhat <- circ.lm$fitted

#Plot the residuals on the unit circle
angle_plot(y - yhat)
circ.lm$coefficients

#In this case, the residuals appear to be symmetrically distributed around 0,
#which is good. We could perform a Kolmogorov-Smirnov test to determine whether
#they are Von Mises distributed. 


#We now check to see if the regression is performed in the way that I think it
#is. To check this, we fit two linear models for cos(y) and sin(y) to the same
#set of predictors.

cos_df <- data.frame(cos(y), cos(x), cos(2*x), sin(x), sin(2*x))
sin_df <- data.frame(sin(y), cos(x), cos(2*x), sin(x), sin(2*x))

colnames(cos_df) <- c("cosy", "cosx", "cos2x", "sinx", "sin2x")
colnames(sin_df) <- c("siny", "cosx", "cos2x", "sinx", "sin2x")
cos.lm <- lm(cosy~.,cos_df)
sin.lm <- lm(siny~.,sin_df)

#Check that coefficients agree
print(cos.lm$coefficients - circ.lm$coefficients[,1])
print(sin.lm$coefficients - circ.lm$coefficients[,2])

#This clearly demonstrates that the circular regression is performed by using
#a linear model for the cosine and sine of the response as a function of the
#cosines and sines of the predictors at different frequencies. It provides a
#clear suggestion for how we should extend this to multiple predictors and a
#clear path to regularization.


#For circular-circular regression, this package does not provide the standard
#error of any of the parameters. We can find them by bootstrapping.
boot_helper<-function(Z,i) {
  Zboot <- Z[i,]
  Y <- Zboot[[1]]
  X <- Zboot[[2]]
  out <- lm.circular(Y,X,type="c-c", order=2)
  params <- out$coefficients
}
simple_boot <- boot(data.frame(y,x), boot_helper, R=2000)
#The entries of enzymes_boot iterate down the left column of circ.lm$coefficients
#and then the right column. So the first entry is the intercept of the cosine
#term, the second is the cos.x1 predictor of the cosine term, etc.

cos_x2_boot <- simple_boot$t[,3]
sin_x2_boot <- simple_boot$t[,5]

#Plot a histogram of the cos(2x) coefficient
ggplot(data.frame(cos_x2_boot), aes(x=cos_x2_boot, y = after_stat(density))) +
  geom_histogram(bins=100) +
  labs(x = "Coefficient",
       y = "Density")

ggplot(data.frame(sin_x2_boot), aes(sin_x2_boot, y = after_stat(density))) +
  geom_histogram(bins=100) +
  labs(x = "Coefficient",
       y = "Density")

#Bootstrapping is not particularly well-behaved in this case! Especially for
#the sine term.


#In this second model, we consider an relationship that is more difficult to
#model using trigonometric polynomials.
x <- circular(runif(200, 0, 2*pi))
y <- 3*x +
  rvonmises(n=50, mu=circular(0), kappa=10)

#Note that x is an odd function, so its complete Fourier series is its Fourier
#sine series. The sine terms should be significant (up to a point) and the
#cosine terms should not be.

circ.lm <- lm.circular(y, x, type="c-c", order=1)
yhat <- circ.lm$fitted
angle_plot(y - yhat)

#With just one term, the model is obviously bad.

circ.lm <- lm.circular(y, x, type="c-c", order=2)
yhat <- circ.lm$fitted
angle_plot(y - yhat)

#Even with two terms, the model is still bad.

circ.lm <- lm.circular(y, x, type="c-c", order=3)
yhat <- circ.lm$fitted
angle_plot(y - yhat)

#With 3 terms, the model suddenly becomes much better, and the residuals
#get much closer to 0. It is clear why we should need exactly 3 terms: the
#y = 3x relationship. If y = 3x then cos(y) = cos(3x) and sin(y) = sin(3x).

circ.lm <- lm.circular(y, x, type="c-c", order=50)
yhat <- circ.lm$fitted
angle_plot(y - yhat)

#With many more terms, the error does not appear to improve that much. It is
#likely that terms beyond 3x are not statistically significant.


#Example 3: infinite number of terms needed
x <- circular(runif(500, 0, 2*pi))
y <- circular(4 + x^2 + rvonmises(n=500, mu=circular(0), kappa=20), modulo="2pi")

circ.lm <- lm.circular(y, x, type="c-c", order=40)
yhat <- circ.lm$fitted
angle_plot(y - yhat)

#For this relationship, where you can always make things more accurate by
#adding more terms, we are at risk of overfitting. 
#We use 10-fold cross-validation to find the order of this model that minimizes
#CV error. This is based on code by Professor Daniel Apley (IEMS).
cv_data <- data.frame(y,x)
num_terms <- seq(1,40,1)
Nrep<-50
K<-10 
n.models = length(num_terms) #number of different models to fit
n=nrow(cv_data)
#y<-nn_df$Y
yhat=matrix(0,n,n.models)
MSE<-matrix(0,Nrep,n.models)
for (j in 1:Nrep) { #Iterate over the number of repetitions
  Ind<-CVInd(n,K)
  for (k in 1:K) { #Iterate over the folds
    for (i in num_terms) { #Iterate over the hyperparameters
      
      #Fit the model
      out <- lm.circular(y[-Ind[[k]]], x[-Ind[[k]]], type="c-c", order=i)
      
      #Predict the hidden fold. We do this manually. This code was written
      #in part by Google Gemini.
      x_new <- x[Ind[[k]]]
      coefs <- out$coefficients
      cos_pred_1 <- lapply(1:i, function (s) coefs[1+s,1] * cos(s*x_new))
      cos_pred_2 <- lapply(1:i, function (s) coefs[1+i+s,1] * sin(s*x_new))
      cos_pred <- Reduce("+", cos_pred_1) + Reduce("+", cos_pred_2) + coefs[1,1]
      sin_pred_1 <- lapply(1:i, function (s) coefs[1+s,2] * cos(s*x_new))
      sin_pred_2 <- lapply(1:i, function (s) coefs[1+i+s,2] * sin(s*x_new))
      sin_pred <- Reduce("+", sin_pred_1) + Reduce("+", sin_pred_2) + coefs[1,2]
      yhat <- circular(atan2(sin_pred, cos_pred), modulo="2pi")
      MSE[j,i] = sum((yhat - y[Ind[[k]]])^2)/n
    }
  }
  print(cat(j, " Repeats done"))
}
MSE
MSEAve<- apply(MSE,2,mean); MSEAve #averaged mean square CV error
MSEsd <- apply(MSE,2,sd); MSEsd   #SD of mean square CV error
r2<-1-MSEAve/var(y); r2  #CV r^2

#The CV results suggest that cross validation MSE is roughly flat once we get
#past 15 terms. We fit the model one last time with 15 cos & sin terms.
circ.lm <- lm.circular(y,x, type="c-c", order=15)
angle_plot(circ.lm$fitted - y)

#Prediction test (make sure this works before CV)
#x <- circular(runif(200, 0, 2*pi))
#y <- 3*x +
#  rvonmises(n=200, mu=circular(0), kappa=10)

#circ.lm <- lm.circular(y, x, type="c-c", order=1)
#yhat <- circ.lm$fitted
#angle_plot(y - yhat)

#coefs <- circ.lm$coefficients

#Order 1 model test
#yhat_cos <- rep(coefs[1,1],length(y)) + coefs[2,1] * cos(x) + coefs[3,1] * sin(x)
#yhat_sin <- rep(coefs[1,2],length(y)) + coefs[2,2] * cos(x) + coefs[3,2] * sin(x)
#yhat_predict <- atan2(yhat_sin, yhat_cos) + 2*pi

#Order 3 model test
#circ.lm <- lm.circular(y, x, type="c-c", order=3)
#yhat <- circ.lm$fitted
#angle_plot(y - yhat)

#coefs <- circ.lm$coefficients

#yhat_cos <- rep(coefs[1,1],length(y)) + coefs[2,1] * cos(x) + coefs[3,1] * cos(2*x) +
#  coefs[4,1] * cos(3*x) + coefs[5,1] * sin(x) + coefs[6,1] * sin(2*x) + coefs[7,1] * sin(3*x)
#yhat_sin <- rep(coefs[1,2],length(y)) + coefs[2,2] * cos(x) + coefs[3,2] * cos(2*x) +
#  coefs[4,2] * cos(3*x) + coefs[5,2] * sin(x) + coefs[6,2] * sin(2*x) + coefs[7,2] * sin(3*x)

#yhat_predict <- atan2(yhat_sin, yhat_cos) + 2*pi #Check


#cos_pred_1 <- lapply(1:i, function (s) coefs[1+s,1] * cos(s*x))
#cos_pred_2 <- lapply(1:i, function (s) coefs[1+i+s,1] * sin(s*x))
#cos_pred <- Reduce("+", cos_pred_1) + Reduce("+", cos_pred_2) + coefs[1,1]

#sin_pred_1 <- lapply(1:i, function (s) coefs[1+s,2] * cos(s*x))
#sin_pred_2 <- lapply(1:i, function (s) coefs[1+i+s,2] * sin(s*x))
#sin_pred <- Reduce("+", sin_pred_1) + Reduce("+", sin_pred_2) + coefs[1,2]
