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

