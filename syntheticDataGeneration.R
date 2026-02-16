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
write.csv(data.frame(y,x), file = "Desktop/NSH_WI26_Rotation" , row.names=FALSE)