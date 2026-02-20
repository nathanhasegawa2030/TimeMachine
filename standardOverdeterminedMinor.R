###standardOverdeterminedMinor.R

#Written by Nathan Hasegawa

#Generates synthetic data where the predictors and responses are uniformly
#distributed, but there are more predictors than responses. Overfitting is a
#serious concern.

#Get helper functions from source files
source("Desktop/NSH_WI26_Rotation/regularized_circular.R")
source("Desktop/NSH_WI26_Rotation/regularized_circular_mgaussian.R")

#Folder name for saving images
folderpath <- "Desktop/NSH_WI26_Rotation/standardOverdeterminedMinorImages"

#Import the training data
train <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/standardOverdeterminedMinorTrain.csv")
y <- train[,1]
x <- train[,-1]

#Plot the distribution of training responses
p <- angle_plot(y, pcol="#2e8b56", title="Distribution of Responses")
p
ggsave(filename="ResponseDistribution.png", plot = p,
       path=folderpath,
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
       path=folderpath,
       width=6,height=6,units="in")
ggsave(filename="MGaussianTrainingResiduals.png", plot = p2,
       path=folderpath,
       width=6,height=6,units="in")
p1
p2

#Import test data
test <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/standardOverdeterminedMinorTest.csv")
y.test <- test[,1]
x.test <- test[,-1]
yhat.test.std <- predict.circ(model.standard, x.test)
yhat.test.mga <- predict.circ.mgaussian(model.mgaussian, x.test)
resid.test.std <- sapply(y.test - yhat.test.std, rescale_angle)
resid.test.mga <- sapply(y.test - yhat.test.mga, rescale_angle)
p1 <- angle_plot(resid.test.std, pcol="#007fff", title="Test Residuals — Standard")
p2 <- angle_plot(resid.test.mga, pcol="#fd1900", title="Test Residuals — MGaussian")
ggsave(filename="StandardTestResiduals.png", plot = p1,
       path=folderpath,
       width=6,height=6,units="in")
ggsave(filename="MGaussianTestResiduals.png", plot = p2,
       path=folderpath,
       width=6,height=6,units="in")
p1
p2


#Show the donut plot for the residuals. This gets at the predicted y-values.
p1 <- donut_plot(yhat.test.std, y.test, title="Test — Standard")
p2 <- donut_plot(yhat.test.mga, y.test, title="Test — MGaussian")
ggsave(filename="StandardTestDonut.png", plot = p1,
       path=folderpath,
       width=6,height=6,units="in")
ggsave(filename="MGaussianTestDonut.png", plot = p2,
       path=folderpath,
       width=6,height=6,units="in")
p1
p2

#Store the fitted models and all the data
save(model.standard, model.mgaussian, test, train, alpha, resid.test.std,
     resid.test.mga, resid.train.std, resid.train.mga,
     file="Desktop/NSH_WI26_Rotation/Workspaces/standardOverdeterminedMinor.RData")

#CAN START HERE

load("Desktop/NSH_WI26_Rotation/Workspaces/standardOverdeterminedMinor.RData")
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

#Filter by where CV r^2 is large, so that we zoom in near the top of the plot

sin.top <- sin.agg[sin.agg$r2 > 0.8,]
cos.top <- cos.agg[cos.agg$r2 > 0.8,]
mga.top <- mga.agg[mga.agg$r2 > 0.8,]

#Plot this again
sin.ref <- sc3D(sin.top)
cos.ref <- sc3D(cos.top)
mga.ref <- sc3D(mga.top)

#Same as above, but now color is nonzero terms and z-axis is CV r^2
sin.all <- sc3D(sin.top, ccol="nonzero")
cos.all <- sc3D(cos.top, ccol="nonzero")
mga.all <- sc3D(mga.top, ccol="nonzero")

#2D version of these plots
sin.ref.2D <- sc2D(sin.top)
cos.ref.2D <- sc2D(cos.top)
mga.ref.2D <- sc2D(mga.top)

#Plot number of nonzero terms in this regime
sin.nzr <- sc3D(sin.top, zcol="nonzero")
cos.nzr <- sc3D(cos.top, zcol="nonzero")
mga.nzr <- sc3D(mga.top, zcol="nonzero")

#2D plot of nonzero terms, confined to the lower regions (less than 10; the
#correct amount is 5)
sin.par <- sin.top[sin.top$nonzero < 10,]
cos.par <- cos.top[cos.top$nonzero < 10,]
mga.par <- mga.top[mga.top$nonzero < 10,]

sin.nzr.2D <- sc2D(sin.par, zcol="nonzero")
cos.nzr.2D <- sc2D(cos.par, zcol="nonzero")
mga.nzr.2D <- sc2D(mga.par, zcol="nonzero")


#Export the plots
folderpath <- "Desktop/NSH_WI26_Rotation/standardOverdeterminedMinorImages"
save_plotly(sin.nzr.2D, folderpath, "SinNonzero.png")
save_plotly(cos.nzr.2D, folderpath, "CosNonzero.png")
save_plotly(mga.nzr.2D, folderpath, "MgaNonzero.png")

save_plotly(sin.ref.2D, folderpath, "Sinr2.png")
save_plotly(cos.ref.2D, folderpath, "Cosr2.png")
save_plotly(mga.ref.2D, folderpath, "Mgar2.png")

