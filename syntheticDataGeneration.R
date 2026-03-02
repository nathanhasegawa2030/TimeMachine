#syntheticDataGeneration.R

#Written by Nathan Hasegawa

#Generates the synthetic training and test data for use in the test cases.
library(circular)
set.seed(33)

#(1) Standard Response Test
xmat <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardTrain.csv", 
          row.names=FALSE)

xmat <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardTest.csv", 
          row.names=FALSE)

#(2) Bilevel Response Test
xmat <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(6*cos(x1) - 14*cos(x2) + 33*sin(x1) + 4*sin(x2),
           0.2*cos(x1) - 0.21*cos(x2) + 0.01*sin(x1) - 0.4*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/bilevelTrain.csv", 
          row.names=FALSE)

xmat <- matrix(runif(25000, 0, 2*pi), nrow = 500, ncol = 50)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(6*cos(x1) - 14*cos(x2) + 33*sin(x1) + 4*sin(x2),
           0.2*cos(x1) - 0.21*cos(x2) + 0.01*sin(x1) - 0.4*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/bilevelTest.csv", 
          row.names=FALSE)

#(3) Nonuniform Predictor Test
xmat <- matrix(rvonmises(25000, (pi/2), 6), nrow = 500, ncol = 50)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(5*cos(x1) - 4*cos(x2) + 2*sin(x1) - 1.9*sin(x2),
           0.5*cos(x1) - 0.75*cos(x2) - 6*sin(x1) + 5.4*sin(x2)) + 
  rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformTrain.csv", 
          row.names=FALSE)

xmat <- matrix(rvonmises(25000, (pi/2), 6), nrow = 500, ncol = 50)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(5*cos(x1) - 4*cos(x2) + 2*sin(x1) - 1.9*sin(x2),
           0.5*cos(x1) - 0.75*cos(x2) - 6*sin(x1) + 5.4*sin(x2)) + 
  rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformTest.csv", 
          row.names=FALSE)

#(4) Standard Overdetermined— Minor
xmat <- matrix(runif(500000, 0, 2*pi), nrow = 500, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardOverdeterminedMinorTrain.csv", 
          row.names=FALSE)

xmat <- matrix(runif(500000, 0, 2*pi), nrow = 500, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=100, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardOverdeterminedMinorTest.csv", 
          row.names=FALSE)

#(5) Standard Overdetermined— Major
xmat <- matrix(runif(100000, 0, 2*pi), nrow = 100, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=100, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardOverdeterminedMajorTrain.csv", 
          row.names=FALSE)

xmat <- matrix(runif(100000, 0, 2*pi), nrow = 100, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=100, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardOverdeterminedMajorTest.csv", 
          row.names=FALSE)

#(6) Nonuniform Overdetermined— Minor
xmat <- matrix(rvonmises(500000, (pi/2), 6), nrow = 500, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(5*cos(x1) - 4*cos(x2) + 2*sin(x1) - 1.9*sin(x2),
           0.5*cos(x1) - 0.75*cos(x2) - 6*sin(x1) + 5.4*sin(x2)) + 
  rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformOverdeterminedMinorTrain.csv", 
          row.names=FALSE)

xmat <- matrix(rvonmises(500000, (pi/2), 6), nrow = 500, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(5*cos(x1) - 4*cos(x2) + 2*sin(x1) - 1.9*sin(x2),
           0.5*cos(x1) - 0.75*cos(x2) - 6*sin(x1) + 5.4*sin(x2)) + 
  rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformOverdeterminedMinorTest.csv", 
          row.names=FALSE)

#(7) Nonuniform Overdetermined— Major
xmat <- matrix(rvonmises(100000, (pi/2), 6), nrow = 100, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(5*cos(x1) - 4*cos(x2) + 2*sin(x1) - 1.9*sin(x2),
           0.5*cos(x1) - 0.75*cos(x2) - 6*sin(x1) + 5.4*sin(x2)) + 
  rvonmises(n=100, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformOverdeterminedMajorTrain.csv", 
          row.names=FALSE)

xmat <- matrix(rvonmises(100000, (pi/2), 6), nrow = 100, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
y <- atan2(5*cos(x1) - 4*cos(x2) + 2*sin(x1) - 1.9*sin(x2),
           0.5*cos(x1) - 0.75*cos(x2) - 6*sin(x1) + 5.4*sin(x2)) + 
  rvonmises(n=100, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/nonuniformOverdeterminedMajorTest.csv", 
          row.names=FALSE)

#(8) Standard Overdetermined— Major, With Correlated Predictors
xmat <- matrix(runif(100000, 0, 2*pi), nrow = 100, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
#Correlated predictors
x[,3] <- (0.75*pi) + x1 + rvonmises(100, mu=0, kappa=250)
x[,4] <- (-0.5*pi) + x2 + rvonmises(100, mu=0, kappa=250)
x[,5] <- (-0.25*pi) + x1 + rvonmises(100, mu=0, kappa=100)
x[,6] <- (pi/6) + x2 + rvonmises(100, mu=0, kappa=100)

y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=100, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardCorrelatedMajorTrain.csv", 
          row.names=FALSE)

xmat <- matrix(runif(100000, 0, 2*pi), nrow = 100, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
#Correlated predictors
x[,3] <- (0.75*pi) + x1 + rvonmises(100, mu=0, kappa=250)
x[,4] <- (-0.5*pi) + x2 + rvonmises(100, mu=0, kappa=250)
x[,5] <- (-0.25*pi) + x1 + rvonmises(100, mu=0, kappa=100)
x[,6] <- (pi/6) + x2 + rvonmises(100, mu=0, kappa=100)
y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=100, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardCorrelatedMajorTest.csv", 
          row.names=FALSE)

#(9) Standard Overdetermined— Minor, With Correlated Predictors
xmat <- matrix(runif(500000, 0, 2*pi), nrow = 500, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
#Correlated predictors
x[,3] <- (0.75*pi) + x1 + rvonmises(500, mu=0, kappa=250)
x[,4] <- (-0.5*pi) + x2 + rvonmises(500, mu=0, kappa=250)
x[,5] <- (-0.25*pi) + x1 + rvonmises(500, mu=0, kappa=100)
x[,6] <- (pi/6) + x2 + rvonmises(500, mu=0, kappa=100)

y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardCorrelatedMinorTrain.csv", 
          row.names=FALSE)

xmat <- matrix(runif(500000, 0, 2*pi), nrow = 500, ncol = 1000)
x <- as.data.frame(xmat)
x1 <- x[,1]
x2 <- x[,2]
#Correlated predictors
x[,3] <- (0.75*pi) + x1 + rvonmises(500, mu=0, kappa=250)
x[,4] <- (-0.5*pi) + x2 + rvonmises(500, mu=0, kappa=250)
x[,5] <- (-0.25*pi) + x1 + rvonmises(500, mu=0, kappa=100)
x[,6] <- (pi/6) + x2 + rvonmises(500, mu=0, kappa=100)

y <- atan2(1.4*cos(x1) - 1.1*cos(x2) + 0.8*sin(x1) + 0.9*sin(x2),
           1.1*cos(x1) - 0.9*cos(x2) - 1.2*sin(x1) - 0.6*sin(x2)) + rvonmises(n=500, mu=circular(0), kappa=100)
write.csv(data.frame(y,x), 
          file = "Desktop/NSH_WI26_Rotation/SyntheticData/standardCorrelatedMinorTest.csv", 
          row.names=FALSE)

