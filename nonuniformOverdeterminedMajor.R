###nonuniformOverdeterminedMajor.R

#Written by Nathan Hasegawa

#Generates synthetic data where the predictors are biased towards the upper
#half of the unit circle, and there are 1000 predictors and 100 training rows.
#Overfitting is a very serious concern here.

#Get helper functions from source files
source("Desktop/NSH_WI26_Rotation/regularized_circular.R")
source("Desktop/NSH_WI26_Rotation/regularized_circular_mgaussian.R")

#Folder name for saving images
folderpath <- "Desktop/NSH_WI26_Rotation/nonuniformOverdeterminedMajorImages"

#Import the training data
train <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformOverdeterminedMajorTrain.csv")
y <- train[,1]
x <- train[,-1]

#Plot the distribution of training responses
p <- angle_plot(y, pcol="#2e8b56", title="Distribution of Responses")
p
ggsave(filename="ResponseDistribution.png", plot = p,
       path=folderpath,
       width=6,height=6,units="in")

#Plot the distribution of x_33 (the x's are all generated according to this
#distribution)
p <- angle_plot(x[,33], pcol="#ffa500", title = "Predictor Distributuion")
p
ggsave(filename="PredictorDistribution.png", plot = p,
       path=folderpath,
       width=6,height=6,units="in")


#Fit both models
alpha = seq(0.5,1,0.05)
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
test <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformOverdeterminedMajorTest.csv")
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
     file="Desktop/NSH_WI26_Rotation/Workspaces/nonuniformOverdeterminedMajor.RData")

#CAN START HERE

load("Desktop/NSH_WI26_Rotation/Workspaces/nonuniformOverdeterminedMajor.RData")
folderpath <- "Desktop/NSH_WI26_Rotation/nonuniformOverdeterminedMajorImages"
source("Desktop/NSH_WI26_Rotation/plotting.R")
source("Desktop/NSH_WI26_Rotation/regularized_circular.R")
source("Desktop/NSH_WI26_Rotation/regularized_circular_mgaussian.R")
test <- read.csv("Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformOverdeterminedMajorTest.csv")

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
sin.top <- sin.agg[sin.agg$r2 > 0.6,]
cos.top <- cos.agg[cos.agg$r2 > 0,]
mga.top <- mga.agg[mga.agg$r2 > 0.5,]

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


#Compute the composite r^2 for the standard case

#Choose the "correct" number of terms that we want to consider. (This is 
#synthetic data; we happen to know the correct number of terms.)
nterms <- c(4,5)
pct <- 0

#Compute the composite r^2's
comp.agg <- comp.r2(model.standard$y, model.standard$x, sin.agg, cos.agg, test,
                    nterms, pct)

#Plot points in blue where the cosine and sine lambdas are equal. These are the
#models that the Mgaussian approach would select from
# minl <- max(min(comp.agg$sinloglambda), min(comp.agg$cosloglambda))
# maxl <- min(max(comp.agg$sinloglambda), max(comp.agg$cosloglambda))
# minz <- min(comp.agg$r2)
# maxz <- max(comp.agg$r2)
# x <- c(minl, minl, maxl, maxl)
# y <- c(minl, minl, maxl, maxl)
# z <- c(minz, maxz, maxz, minz)
# bluept <- data.frame(x,y,z)

#3D plot of CV r^2 for the composite model

comp.r2.p <- sc3D(comp.agg, xcol="sinloglambda", ycol="cosloglambda",
                  xlab='log<sub>10</sub>(λ<sub>sin</sub>)',
                  ylab='log<sub>10</sub>(λ<sub>cos</sub>)')


#Now plot the test r^2 too
test.r2.p <- sc3D(comp.agg, xcol="sinloglambda", ycol="cosloglambda",
                  zcol="testr2",
                  xlab='log<sub>10</sub>(λ<sub>sin</sub>)',
                  ylab='log<sub>10</sub>(λ<sub>cos</sub>)',
                  zlab='Test r<sup>2</sup>')

#Now do the same for the overfit models
nterms <- seq(1,2*dim(model.standard$x)[2],1)
pct <- 0.95
comp.overfit <- comp.r2(model.standard$y, model.standard$x, sin.agg, cos.agg, test,
                        nterms, pct)

#3D plot of CV r^2 for the overfit models
comp.overfit.r2.p <- sc3D(comp.overfit, xcol="sinloglambda", ycol="cosloglambda",
                          xlab='log<sub>10</sub>(λ<sub>sin</sub>)',
                          ylab='log<sub>10</sub>(λ<sub>cos</sub>)')

#Same for test r^2
test.overfit.r2.p <- sc3D(comp.overfit, xcol="sinloglambda", ycol="cosloglambda",
                          zcol="testr2",
                          xlab='log<sub>10</sub>(λ<sub>sin</sub>)',
                          ylab='log<sub>10</sub>(λ<sub>cos</sub>)',
                          zlab='Test r<sup>2</sup>')

#2D versions of the above plots

comp.r2.2D <- sc2D(comp.agg, xcol="cosloglambda", ycol="sinloglambda",
                   xlab='log<sub>10</sub>(λ<sub>cos</sub>)',
                   ylab='log<sub>10</sub>(λ<sub>sin</sub>)')

test.r2.2D <- sc2D(comp.agg, xcol="cosloglambda", ycol="sinloglambda",
                   zcol="testr2",
                   xlab='log<sub>10</sub>(λ<sub>cos</sub>)',
                   ylab='log<sub>10</sub>(λ<sub>sin</sub>)',
                   zlab='Test r<sup>2</sup>')

comp.overfit.2D <- sc2D(comp.overfit, xcol="cosloglambda", ycol="sinloglambda",
                        xlab='log<sub>10</sub>(λ<sub>cos</sub>)',
                        ylab='log<sub>10</sub>(λ<sub>sin</sub>)')

test.overfit.2D <- sc2D(comp.overfit, xcol="cosloglambda", ycol="sinloglambda",
                        zcol="testr2",
                        xlab='log<sub>10</sub>(λ<sub>cos</sub>)',
                        ylab='log<sub>10</sub>(λ<sub>sin</sub>)',
                        zlab='Test r<sup>2</sup>')

#Plot of lambda vs. composite r^2 for MGaussian. Here we use the entire range of
#lambda's, since there aren't that many, and prediction is fast.
mga.a1 <- aggregate(mga.data[mga.data$alpha == 1,])
mga.a1.plot <- mga.a1[mga.a1$r2 > 0.5,]
mga.r2.2D <- sc2D(mga.a1.plot, xcol="loglambda", ycol="r2", zcol="nonzero",
                  ylab='CV r<sup>2</sup>')
#Same plot but for sparse models only
mga.a1.plot <- mga.a1.plot[mga.a1.plot$nonzero < 10,]
mga.r2.2D.sparse <- sc2D(mga.a1.plot, xcol="loglambda", ycol="r2", zcol="nonzero",
                  ylab='CV r<sup>2</sup>')


#Plot of lambda vs. test r^2 for MGaussian
mga.model <- model.mgaussian$model
mga.test.pred <- predict.mgaussian(mga.model, test[,-1])
y.test <- rescale_angle(test[,1])
resid2 <- function (ypred) {mean((rescale_angle(ypred - y.test))^2)}
test.resid.mean <- sapply(mga.test.pred, resid2)
mga.test.r2 <- 1 - test.resid.mean/(var(cos(y.test))+var(sin(y.test)))
mga.a1$testr2 <- mga.test.r2
mga.a1.plot <- mga.a1[mga.a1$testr2 > 0.4,]
test.mga.2D <- sc2D(mga.a1.plot, xcol="loglambda", ycol="testr2", zcol="nonzero",
                    ylab='Test r<sup>2</sup>')
#Same plot but for sparse models only
mga.a1.plot <- mga.a1.plot[mga.a1.plot$nonzero < 10,]
test.mga.2D.sparse <- sc2D(mga.a1.plot, xcol="loglambda", ycol="testr2", zcol="nonzero",
                         ylab='Test r<sup>2</sup>')


#Export the plots
save_plotly(sin.nzr.2D, folderpath, "SinNonzero.png")
save_plotly(cos.nzr.2D, folderpath, "CosNonzero.png")
save_plotly(mga.nzr.2D, folderpath, "MgaNonzero.png")

save_plotly(sin.ref.2D, folderpath, "Sinr2.png")
save_plotly(cos.ref.2D, folderpath, "Cosr2.png")
save_plotly(mga.ref.2D, folderpath, "Mgar2.png")

save_plotly(comp.r2.2D, folderpath, "CompositeSparser2.png")
save_plotly(comp.overfit.2D, folderpath, "CompositeOverfitr2.png")
save_plotly(test.r2.2D, folderpath, "CompositeTestSparser2.png")
save_plotly(test.overfit.2D, folderpath, "CompositeTestOverfitr2.png")
save_plotly(mga.r2.2D, folderpath, "Mgar2.png")
save_plotly(test.mga.2D, folderpath, "MgaTestr2.png")
save_plotly(mga.r2.2D.sparse, folderpath, "Mgar2Sparse.png")
save_plotly(test.mga.2D.sparse, folderpath, "MgaTestr2Sparse.png")
