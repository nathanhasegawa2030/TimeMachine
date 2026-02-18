###nonuniformResponseTest.R

#Written by Nathan Hasegawa

#Generates synthetic data where the predictors are NOT uniformly distributed but
#the responses are well spread out across the unit circle.

#Get helper functions from source files
source("Desktop/NSH_WI26_Rotation/regularized_circular.R")
source("Desktop/NSH_WI26_Rotation/regularized_circular_mgaussian.R")

#Import the training data
train <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformTrain.csv")
y <- train[,1]
x <- train[,-1]

#Plot the distribution of x_33 (the x's are all generated according to this
#distribution)
p <- angle_plot(x[,33], pcol="#ffa500", title = "Predictor Distributuion")
p
ggsave(filename="PredictorDistribution.png", plot = p,
       path="Desktop/NSH_WI26_Rotation/nonuniformResponseTestImages",
       width=6,height=6,units="in")

#Plot the distribution of training responses
p <- angle_plot(y, pcol="#2e8b56", title="Distribution of Responses")
p
ggsave(filename="ResponseDistribution.png", plot = p,
       path="Desktop/NSH_WI26_Rotation/nonuniformResponseTestImages",
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
       path="Desktop/NSH_WI26_Rotation/nonuniformResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTrainingResiduals.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/nonuniformResponseTestImages",
       width=6,height=6,units="in")
p1
p2

#Import test data
test <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformTest.csv")
y.test <- test[,1]
x.test <- test[,-1]
yhat.test.std <- predict.circ(model.standard, x.test)
yhat.test.mga <- predict.circ.mgaussian(model.mgaussian, x.test)
resid.test.std <- sapply(y.test - yhat.test.std, rescale_angle)
resid.test.mga <- sapply(y.test - yhat.test.mga, rescale_angle)
p1 <- angle_plot(resid.test.std, pcol="#007fff", title="Test Residuals — Standard")
p2 <- angle_plot(resid.test.mga, pcol="#fd1900", title="Test Residuals — MGaussian")
ggsave(filename="StandardTestResiduals.png", plot = p1,
       path="Desktop/NSH_WI26_Rotation/nonuniformResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTestResiduals.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/nonuniformResponseTestImages",
       width=6,height=6,units="in")
p1
p2


#Show the donut plot for the residuals. This gets at the predicted y-values.
p1 <- donut_plot(yhat.test.std, y.test, title="Test — Standard")
p2 <- donut_plot(yhat.test.mga, y.test, title="Test — MGaussian")
ggsave(filename="StandardTestDonut.png", plot = p1,
       path="Desktop/NSH_WI26_Rotation/nonuniformResponseTestImages",
       width=6,height=6,units="in")
ggsave(filename="MGaussianTestDonut.png", plot = p2,
       path="Desktop/NSH_WI26_Rotation/nonuniformResponseTestImages",
       width=6,height=6,units="in")
p1
p2

#Store the fitted models and all the data
save(model.standard, model.mgaussian, test, train, alpha, resid.test.std,
     resid.test.mga, resid.train.std, resid.train.mga, aggregate, angle_plot,
     donut_plot, rescale_angle,
     file="Desktop/NSH_WI26_Rotation/Workspaces/nonuniformResponseTest.RData")
load("Desktop/NSH_WI26_Rotation/Workspaces/nonuniformResponseTest.RData")

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


#3D Plots of CV r^2 vs. alpha and lambda. We exclude alpha = 0 from this plot.
sin.p <- plot_ly(sin.agg, x = ~loglambda, y = ~alpha, z = ~r2)
sin.p <- sin.p %>% layout(scene = list(xaxis = list(title = 'Log10 Lambda'),
                                       yaxis = list(title = 'Alpha', range = c(0.05, 1)),
                                       zaxis = list(title = 'CV r^2'), range = c(0.8,1)),
                          font = list(family = "Trebuchet MS", size = 14, 
                                      color = "black"))

cos.p <- plot_ly(cos.agg, x = ~loglambda, y = ~alpha, z = ~r2)
cos.p <- cos.p %>% layout(scene = list(xaxis = list(title = 'Log10 Lambda'),
                                       yaxis = list(title = 'Alpha', range = c(0.05, 1)),
                                       zaxis = list(title = 'CV r^2'), range = c(0.8,1)),
                          font = list(family = "Trebuchet MS", size = 14, 
                                      color = "black"))

mga.p <- plot_ly(mga.agg, x = ~loglambda, y = ~alpha, z = ~r2)
mga.p <- mga.p %>% layout(scene = list(xaxis = list(title = 'Log10 Lambda'),
                                       yaxis = list(title = 'Alpha', range = c(0.05, 1)),
                                       zaxis = list(title = 'CV r^2'), range = c(0.8,1)),
                          font = list(family = "Trebuchet MS", size = 14, 
                                      color = "black"))



ggplot(sin.agg, aes(x=loglambda, y=alpha, color=r2)) + 
  geom_point() +
  scale_color_gradient(low = "#007fff", high = "#fd1900") +
  xlim(min(sin.agg$loglambda)-0.1,0.1)