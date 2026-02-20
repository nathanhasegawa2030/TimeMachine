###standardResponseTest.R

#Written by Nathan Hasegawa

#Imports synthetic data where the predictors are uniformly distributed and the
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
save(model.standard, model.mgaussian, test, train, alpha, resid.test.std,
     resid.test.mga, resid.train.std, resid.train.mga,
     file="Desktop/NSH_WI26_Rotation/Workspaces/standardResponseTest.RData")

#CAN START HERE

load("Desktop/NSH_WI26_Rotation/Workspaces/standardResponseTest.RData")
source("Desktop/NSH_WI26_Rotation/plotting.R")

#Summary statistics for test residuals
summary(resid.test.std)
summary(resid.test.mga)

#Summary statistics for training residuals
summary(resid.train.std)
summary(resid.train.mga)

#Plot of CV error as function of alpha and lambda. Google Gemini helped
#write this code
sin.data <- model.standard$sin.data
cos.data <- model.standard$cos.data
mga.data <- model.mgaussian$cv.data

#Take log10 of lambda for interpolation
sin.data$loglambda = log10(sin.data$lambda)
cos.data$loglambda = log10(cos.data$lambda)
mga.data$loglambda = log10(mga.data$lambda)

#Take the average of responses with equal x- and y-components
sin.agg <- aggregate(sin.data)
cos.agg <- aggregate(cos.data)
mga.agg <- aggregate(mga.data)

#Filter by where alpha > 0. We do this because when alpha = 0, there is no
#LASSO penalty, and terms do not go to 0. That is a bad model.
sin.agg <- sin.agg[sin.agg$alpha > 0,]
cos.agg <- cos.agg[cos.agg$alpha > 0,]
mga.agg <- mga.agg[mga.agg$alpha > 0,]

#3D Plots of CV r^2 vs. alpha and lambda
sin.p <- sc3D(sin.agg)
cos.p <- sc3D(cos.agg)
mga.p <- sc3D(mga.agg)

#Filter by where CV r^2 is large for best results. TO DEFINE "LARGE," WE HAVE TO
#LOOK AT THE PLOT FIRST.

#For this test case, 0.8 is a good benchmark
sin.top <- sin.agg[sin.agg$r2 > 0.8,]
cos.top <- cos.agg[cos.agg$r2 > 0.8,]
mga.top <- mga.agg[mga.agg$r2 > 0.8,]
  
#Plot this again
sin.ref <- sc3D(sin.top)
cos.ref <- sc3D(cos.top)
mga.ref <- sc3D(mga.top)

#Plot number of nonzero terms in this regime
sin.nzr <- sc3D(sin.top, zcol="nonzero")
cos.nzr <- sc3D(cos.top, zcol="nonzero")
mga.nzr <- sc3D(mga.top, zcol="nonzero")



