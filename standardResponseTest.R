###standardResponseTest.R

#Written by Nathan Hasegawa

#Generates synthetic data where the predictors are uniformly distributed and the
#responses tend to be roughly uniformly distributed. This is a
#test of which method performs better when both responses have very different 
#orders of magnitude.

#Get helper functions from source files
source("Desktop/NSH_WI26_Rotation/regularized_circular.R")
source("Desktop/NSH_WI26_Rotation/regularized_circular_mgaussian.R")

#Import the training data
train <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/standardTrain.csv")
y <- train[,1]
x <- train[,-1]

#Plot the distribution of training responses
p <- angle_plot(y, pcol="#2e8b56", title="Distribution of Responses")
p
ggsave(filename="ResponseDistribution.png", plot = p,
       path="Desktop/NSH_WI26_Rotation/standardResponseTestImages",
       width=6,height=6,units="in")

#Fit both models
alpha = seq(0,1,0.1)
model.standard <- circ.lm(y,x,a=alpha)
model.mgaussian <- circ.lm.mgaussian(y,x,a=alpha)

resid.train.std = model.standard$residuals
resid.train.mga = model.mgaussian$residuals

#Plot the residuals
p1 <- angle_plot(resid.train.std, pcol="#007fff", 
                 title="Training Residuals — Standard")
p2 <- angle_plot(resid.train.mga, pcol="#fd1900", 
                 title="Training Residuals — MGaussian")
ggsave(filename="StandardTrainingResiduals.png", plot = p1,
       path="Desktop/NSH_WI26_Rotation/standardResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTrainingResiduals.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/standardResponseTestImages",
       width=6,height=6,units="in")
p1
p2

#Import test data
test <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/standardTest.csv")
y.test <- test[,1]
x.test <- test[,-1]
yhat.test.std <- predict.circ(model.standard, x.test)
yhat.test.mga <- predict.circ.mgaussian(model.mgaussian, x.test)
resid.test.std <- sapply(y.test - yhat.test.std, rescale_angle)
resid.test.mga <- sapply(y.test - yhat.test.mga, rescale_angle)
p1 <- angle_plot(resid.test.std, pcol="#007fff", title="Test Residuals — Standard")
p2 <- angle_plot(resid.test.mga, pcol="#fd1900", title="Test Residuals — MGaussian")
ggsave(filename="StandardTestResiduals.png", plot = p1,
       path="Desktop/NSH_WI26_Rotation/standardResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTestResiduals.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/standardResponseTestImages",
       width=6,height=6,units="in")
p1
p2


#Show the donut plot for the residuals. This gets at the predicted y-values.
p1 <- donut_plot(yhat.test.std, y.test, title="Test — Standard")
p2 <- donut_plot(yhat.test.mga, y.test, title="Test — MGaussian")
ggsave(filename="StandardTestDonut.png", plot = p1,
       path="Desktop/NSH_WI26_Rotation/standardResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTestDonut.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/standardResponseTestImages",
       width=6,height=6,units="in")
p1
p2

#Store the fitted models and all the data

#Summary statistics for test residuals
summary(resid.test.std)
summary(resid.test.mga)

#Summary statistics for training residuals
summary(model.standard$residuals)
summary(model.mgaussian$residuals)



