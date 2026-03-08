#TMEDA.R

#Authored by Nathan Hasegawa

#Exploratory data analysis for the TimeMachine data. Mainly trying to understand
#the underlying distribution of predictors and responses, and whether the
#distribution may lend itself well to c-c regression.

library(patchwork)

#Load the data and source files. We use only the Z-Score data for this.
source("Desktop/NSH_WI26_Rotation/plotting.R")
load("Desktop/NSH_WI26_Rotation/data/TMexampleData.Rdata")
load("Desktop/NSH_WI26_Rotation/cycling_list.RData")

#Get the Z-score expression data. This code is copied directly from Figure1.R in
#the original TimeMachine repository.

all.expr[is.na(all.expr)] <- 0
cycling.expr <- all.expr[rownames(all.expr) %in% cycling.list,]
sample.mean <- apply(cycling.expr,2,mean)
sample.std <- apply(cycling.expr,2,sd)
predictor.expr.nor <- apply(cycling.expr, 1, function(x) {
  tmp <- (x-sample.mean)/sample.std
})
predictor.expr.nor <- t(predictor.expr.nor) #These are values of predictor variables
#in each sample. NOT the predictor itself (the 37 genes)
predictor.expr.final <- predictor.expr.nor

#There are 50 samples for which we do not have ANY of the predictors (if we are
#using melatonin)
#which(is.na(predictor.expr.final[1,])) == which(is.na(predictor.expr.final[3,]))

#Convert to [-pi,pi] for angle plot
rescale_pi <- function (x) {
  rescale_angle((x*pi) / max(abs(x)))
}
predictor.expr.circular <- t(apply(predictor.expr.final,2,rescale_pi))

#Check that the maximum absolute value in each patient is pi
for (i in 1:1115) {
  if (abs(max(abs(predictor.expr.circular[i,])) - pi) > 0.00001) {
    print("WARNING: Standardization is incorrect")
  }
}

#Plot the distribution of the predictors. Code written with help from Google
#Gemini. Since 37 is prime, we only plot the first 36.
genenames <- rownames(predictor.expr.circular)
predictor_plots <- list()
for (i in 1:36) {
  p = angle_plot(predictor.expr.circular[i,], title=genenames[i])
  predictor_plots[[i]] = p
}
p <- wrap_plots(predictor_plots, nrow = 4, ncol = 9)
ggsave(filename="PredictorDistribution.png", plot = p,
       path="Desktop/NSH_WI26_Rotation/TMEDAImages",
       width=24,height=13,units="in")

#Plot the distribution of the responses. For melatonin, we omit rows where DLMO
#is not available.

DLMO_angles = all.meta$DLMO25 * (pi/12)
DLMO_angles <- DLMO_angles[!is.na(DLMO_angles)]
response_DLMO <- rescale_angle(DLMO_angles)
p <- angle_plot(response_DLMO, title="Response DLMO25", pcol="#2e8b56", mode="clock")
ggsave(filename="DLMODistribution.png", plot = p,
       path="Desktop/NSH_WI26_Rotation/TMEDAImages",
       width=6,height=6,units="in")

response_local_time = rescale_angle(all.meta$LocalTime * (pi/12))
p <- angle_plot(response_local_time, title="Local Time", pcol="#2e8b56", mode="clock")
ggsave(filename="LocalTimeDistribution.png", plot = p,
       path="Desktop/NSH_WI26_Rotation/TMEDAImages",
       width=6,height=6,units="in")


#Do this again for training data ONLY



