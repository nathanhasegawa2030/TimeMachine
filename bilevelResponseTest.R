###bilevelResponseTest.R

#Written by Nathan Hasegawa

#Generates synthetic data where the predictors are uniformly distributed and the
#responses tend to be concentrated towards angles of ~pi/2 and ~3pi/2. This is a
#test of which method performs better when both responses have very different 
#orders of magnitude.

#Get helper functions from source files
source("Desktop/NSH_WI26_Rotation/regularized_circular.R")
source("Desktop/NSH_WI26_Rotation/regularized_circular_mgaussian.R")

#Generate the training data. Note that we have 50 predictors, and only two of
#them are actually used to generate the synthetic responses; this is also a test
#of how good the regularization is.

xmat <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(6*cos(x1) - 14*cos(x2) + 33*sin(x1) + 4*sin(x2),
           0.2*cos(x1) - 0.21*cos(x2) + 0.01*sin(x1) - 0.4*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)

#The plot below shows the distribution of training responses.
#Notice how it is concentrated towards the top and bottom. We can treat this in
#some respects as a classification problem.
p <- angle_plot(y, pcol="#2e8b56", title="Distribution of Responses")
ggsave(filename="ResponseDistribution.png", plot = p,
       path="Desktop/NSH_WI26_Rotation/bilevelResponseTestImages",
       width=6,height=6,units="in")

#Fit both models
alpha = seq(0,1,0.1)

model.standard <- circ.lm(y,x,a=alpha)
model.mgaussian <- circ.lm.mgaussian(y,x,a=alpha)

#Plot the residuals
p1 <- angle_plot(model.standard$residuals, pcol="#007fff", 
                 title="Training Residuals — Standard")
p2 <- angle_plot(model.standard$residuals, pcol="#fd1900", 
                 title="Training Residuals — MGaussian")
ggsave(filename="StandardTrainingResiduals.png", plot = p1,
       path="Desktop/NSH_WI26_Rotation/bilevelResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTrainingResiduals.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/bilevelResponseTestImages",
       width=6,height=6,units="in")

#Due to the bimodal nature of the response distribution, we can treat this in
#essence as a classification problem. We will define the classification accuracy
#as the proportion of responses for which the residual is less than pi/4.
#We compute the training classification accuracy for both models.
training.accuracy.std <- sum(abs(model.standard$residuals) < (pi/4))/length(y)
training.accuracy.mga <- sum(abs(model.mgaussian$residuals) < (pi/4))/length(y)


#We now generate test data to evaluate how good our fitted models are
#at predicting new responses.

#Generate test data
xmat.test <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x.test <- as.data.frame(xmat.test)
x1.t <- x.test[,1]
x2.t <- x.test[,2]
y.test <- atan2(6*cos(x1.t) - 14*cos(x2.t) + 33*sin(x1.t) + 4*sin(x2.t),
                0.2*cos(x1.t) - 0.21*cos(x2.t) + 0.01*sin(x1.t) - 0.4*sin(x2.t)) + 
          rvonmises(n=500, mu=circular(0), kappa=100)

yhat.test.std <- predict.circ(model.standard, x.test)
yhat.test.mga <- predict.circ.mgaussian(model.mgaussian, x.test)
resid.test.std <- sapply(y.test - yhat.test.std, rescale_angle)
resid.test.mga <- sapply(y.test - yhat.test.mga, rescale_angle)
p1 <- angle_plot(resid.test.std, pcol="#007fff", title="Test Residuals — Standard")
p2 <- angle_plot(resid.test.mga, pcol="#fd1900", title="Test Residuals — MGaussian")
ggsave(filename="StandardTestResiduals.png", plot = p1,
       path="Desktop/NSH_WI26_Rotation/bilevelResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTestResiduals.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/bilevelResponseTestImages",
       width=6,height=6,units="in")


#Show the donut plot for the residuals. This gets at the predicted y-values.
p1 <- donut_plot(yhat.test.std, y.test, title="Test — Standard")
p2 <- donut_plot(yhat.test.mga, y.test, title="Test — MGaussian")
ggsave(filename="StandardTestDonut.png", plot = p1,
       path="Desktop/NSH_WI26_Rotation/bilevelResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTestDonut.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/bilevelResponseTestImages",
       width=6,height=6,units="in")

#Compute the classification accuracy for the test data
test.accuracy.std <- sum(abs(resid.test.std) < (pi/4))/length(y)
test.accuracy.mga <- sum(abs(resid.test.mga) < (pi/4))/length(y)
