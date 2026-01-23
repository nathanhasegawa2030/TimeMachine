###Citation for circglmbayes package: https://github.com/keesmulder/circglmbayes/tree/master
###https://www.sciencedirect.com/science/article/abs/pii/S0022249616300827
###Kees Mulder

###This file is the TimeStampFns.R file, modified to perform univariate
###circular regression instead of the XY-plane regression done for TimeMachine.

CXX_STD = C++11

#=====================================================================
# time conversion functions
#=====================================================================

time2dectime <- function(time){
	if(any(c("POSIXct","POSIXt")%in%class(time))){
		# if we have POSIX input, make it a string to meet our (bad)
		# assumptions later...
		time <- strftime(time,"%H:%M")
	}
	if (is.numeric(time)) {
		# horrible unchecked assumption that if input is numeric 
		# it's already decimal 24h time and can be returned
		dectime <- time
	} else {
		# horrible unchecked assumption that otherwise we have
		# times in %H:%M string format -- note that 2:30PM will fail!
		hr_min <- strsplit(as.character(time),":")
		hr <- as.numeric(sapply(hr_min,`[`,1))
		min <- as.numeric(sapply(hr_min,`[`,2))
        min[is.na(min)] <- 0
        dectime <- hr+min/60
	}
	return(dectime)
}


time2angle <- function(time){
	dectime <- time2dectime(time)
	timeAngle <- (dectime%%24)*2*pi/24
	return(timeAngle)
}


time2XY <- function(time){
	a <- time2angle(time)
	XYtime <- zapsmall(cbind(
		timeX=sin(a),
		timeY=cos(a)
	))
	return(XYtime)
}

XY2dectime <- function(XYtime){
	(atan2(XYtime[,1],XYtime[,2])%%(2*pi))*(24/(2*pi))
}

XY2amp <- function(XYtime){
  (sqrt(XYtime[,1]**2 + XYtime[,2]**2))
}

#=====================================================================
# Splitting & centering the expression data for each subject
#=====================================================================

recalibrateExprs <- function(exprMat,subjectIDs){
	# exprMat is a genes * subjects matrix
	# subjectIDs a vector of subject IDs
	exprs <- as.data.frame(t(exprMat))
	exprs <- split(exprs,subjectIDs)
	FCs <- lapply(exprs,scale,center=T,scale=F)
	FCs <- t(do.call(rbind,FCs))[,colnames(exprMat)]
	return(FCs)
}	

#=====================================================================
# FNS FOR CALCULATING PREDICTION ERROR MODULO 24
#=====================================================================

timeErr <- function(trueTimes,predTimes,...){ 
    predTimes <- predTimes%%24
    trueTimes <- trueTimes%%24
    timeDiffs <- abs(predTimes-trueTimes)
    timeDiffs <- pmin(timeDiffs,24-timeDiffs)
    return(timeDiffs)
}

timeErrSigned <- function(trueTimes,predTimes,...){
	predTimes <- predTimes%%24
	trueTimes <- trueTimes%%24
	timeDiffs <- predTimes-trueTimes
	timeDiffs[which(timeDiffs>12)]<-timeDiffs[which(timeDiffs>12)]-24
	timeDiffs[which(timeDiffs<(-12))]<-timeDiffs[which(timeDiffs<(-12))]+24	
	return(timeDiffs)
}


#=====================================================================
# FNS FOR TRAINING AND APPLYING THE CIRCULAR REGRESSION PREDICTOR
#=====================================================================

#Note: no regularization for now. Regularization is needed to ensure
#that coefficients are in (-pi, pi].
stopifnot(require(circular))
stopifnot(require(circglmbayes))

#Compute the model prediction for time. This uses the circular package.
x <- cbind(rnorm(10), rep(1, 10))
x <- cbind(rnorm(10), rep(1,10))
y <- circular(2*atan(c(x%*%c(5,1))), modulo="2pi")+rvonmises(10, mu=circular(0), kappa=100)
model <- lm.circular(y=y, x=x, type="c-l", init=c(5,1), verbose=TRUE)

#Compute the 

expr = CPtrainZScoreDat
subjIDs = CPtrainSubjs
times = CPtrainTimes
trainFrac = 1

x <- t(expr)
y <- time2angle(times) - pi

#NAIVE objective function: response is modded by 2*pi, linear link function.
obj_fcn_naive <- function(beta) {
    #The last element of beta is the constant term; each one before
    #corresponds to one predictor.
    sq_error <- 0
    for (i in seq(1,dim(x)[1])) { #Iterates over rows of x
        #Compute the predictor
        y_pred <- beta[1:dim(x)[2]] %*% x[i,] + beta[length(beta)]
        y_pred <- y_pred %% (2*pi)

        #Compute the squared residual
        resid <- min(abs(y_pred - y[i]), 2*pi - abs(y_pred - y[i]))
        sq_error <- sq_error + resid^2
    }
    #TEMPORARY: Ridge penalty
    reg_penalty <- 0.1 * sum(beta^2)
    return(sq_error + reg_penalty)
}

#Conduct minimization using nlm()
out_nlm_naive <- nlminb(rep(0, dim(x)[2]+1), obj_fcn_naive,
           lower=rep(-pi, dim(x)[2]+1), 
           upper=rep(pi, dim(x)[2]+1))

#This method is very sensitive to the initial condition and isn't great.
#The SSE is not that much below setting all the parameters to 0.

#We now use the link function proposed by Fisher & Lee (1992) and see if that
#works better. The distributional parameter lambda needs to be tuned via 
#cross-validation and is fixed for the purpose of this test.
obj_fcn_tanlink <- function(beta) {
    #The last element of beta is the constant term; each one before
    #corresponds to one predictor.
    sq_error <- 0
    for (i in seq(1,dim(x)[1])) { #Iterates over rows of x
        #Compute the predictor
        g_comp <- beta[1:dim(x)[2]] %*% x[i,] + beta[length(beta)]
        y_pred <- 2 * atan(sign(g_comp) * abs(g_comp)^0.2)

        #Compute the squared residual
        resid <- min(abs(y_pred - y[i]), 2*pi - abs(y_pred - y[i]))
        sq_error <- sq_error + resid^2
    }
    #TEMPORARY: Ridge penalty with lambda = 0.1
    reg_penalty <- 1 * sum(abs(beta))
    return(sq_error + reg_penalty)
}

#Note that with the link function, there is no longer a good reason to bound
#the coefficients by +/- pi
out_nlm_tanlink <- nlm(obj_fcn_tanlink, p=rep(0, dim(x)[2]+1))


#Using circglmbayes package




#Using circular package
test_model <- lm.circular(y=y,x=x,type="c-l",init=rep(-0.01,dim(x)[2]), verbose=TRUE)


#TM.CPhrs.Ratio <- trainTimeStamp(
#  expr=CPtrainRatioDat, 
#  subjIDs=CPtrainSubjs,
#  times=CPtrainTimes,
#  trainFrac=1, 
#  recalib=FALSE,
#  a = 0.2, s = 0.31,
#  foldid=CPtrain.foldid, 
#  plot=FALSE 
#)

#fn <- function(p) {Yhat<-(p[1]*X)/(p[2]+X); sum((Y-Yhat)^2)} 
#out_nlm<-nlm(fn,p=c(gamma_0_guess, gamma_1_guess),hessian=TRUE)
#print(out_nlm)

#=====================================================================
# FNS FOR MODELING X&Y CLOCK COORDS
#=====================================================================

trainTimeStamp <- function(expr,subjIDs,times,trainFrac=0.5,a=NULL,s=NULL,plot=FALSE,recalib=FALSE,...){
	stopifnot(require(glmnet))
	out <- list()

	subjects <- unique(subjIDs)
	trainSubj <- sample(subjects,size=round(trainFrac*length(subjects)),replace=FALSE)
	train <- subjIDs%in%trainSubj
	out$train <- train
	x <- t(expr)
	y <- time2angle(times) #Note: in base TimeMachine, this is time2XY. We are using only
    #the angle (for now).

    #For now, we stick to elastic net.

	out$cv.fit <- cv.glmnet(x[train,],y[train,],keep=T,alpha=a,family="mgaussian",...)
	out$coef <- as.matrix(do.call(cbind,coef(out$cv.fit, s=s)))
	out$coef <- out$coef[rowSums(out$coef!=0)>0,][-1,]
	out$pred <- XY2dectime(predict(out$cv.fit,x,s=s)[,,1])
	if(plot){
		dectime <- time2dectime(times)
		errplot(dectime,out$pred,col=2-train)
	}
	return(out)
}

predTimeStamp <- function(timestamp,newx=NULL,s=NULL){
	stopifnot(require(glmnet))
	if(is.null(newx)){
		return(timestamp$pred)
	}else{
		newx <- t(newx)
		XY2dectime(predict(timestamp$cv.fit,newx,s=s)[,,1])
	}
}

predTimeStampAmp <- function(timestamp,newx=NULL,s=NULL){
  stopifnot(require(glmnet))
  if(is.null(newx)){
    return(timestamp$pred)
  }else{
    newx <- t(newx)
    XY2amp(predict(timestamp$cv.fit,newx,s=s)[,,1])
  }
}

#=====================================================================
# FNS FOR PLOTTING PREDICTED CLOCK COORDS
#=====================================================================

errplot <- function(trueHr,predHr,col=1,pch=1,main="Time of Day (24h)"){
	# plot predicted vs. true time of day , with grey bands at 2 & 4h MOE
	plot(trueHr,predHr,xlab="",ylab="",main=main,yaxs='i', xaxp=c(0,24,6),yaxp=c(0,24,6),col=col,pch=pch,xlim=c(0,24),ylim=c(0,24))
	trueHr <- (trueHr+240)%%24
	predHr <- (predHr+240)%%24
	mtext("True",side=1,line=2.2,cex=0.8)
	mtext("Predicted",side=2,line=2.4,cex=0.8)
	abline(a=0,b=1,col="grey");
	abline(a=24,b=1,col="grey");
	abline(a=-24,b=1,col="grey");

	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-4,-4,4,4),border="grey",lty=3,col=adjustcolor("grey",alpha=0.2));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-2,-2,2,2),border="grey",lty=2,col=adjustcolor("grey",alpha=0.3));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-4,-4,4,4)+24,border="grey",lty=3,col=adjustcolor("grey",alpha=0.2));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-2,-2,2,2)+24,border="grey",lty=2,col=adjustcolor("grey",alpha=0.3));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-4,-4,4,4)-24,border="grey",lty=3,col=adjustcolor("grey",alpha=0.2));
	polygon(c(-48,48,48,-48),c(-48,48,48,-48)+c(-2,-2,2,2)-24,border="grey",lty=2,col=adjustcolor("grey",alpha=0.3));
	points(trueHr,predHr,col=col,pch=pch)
}

tolplot <- function(trueHr,predHr,add=FALSE,col=1){
	# plot % correct by tolerance, dropping lines at median, 90quantile, & 1sd 
	predTimes <- predHr%%24
	trueTimes <- trueHr%%24
	hrerr <- abs(predTimes-trueTimes)
	hrerr <- hrerr[!is.na(hrerr)]
	hrerr <- pmin(hrerr,24-hrerr)
	hrsoff <- seq(0,12,length=49)
	fracacc <- sapply(hrsoff,function(hrtol){
		100*sum(abs(hrerr)>hrtol)/length(hrerr)
	})
	if(!add){
		col=1
		plot(hrsoff,100-fracacc,xlim=c(0,12),
			main="Absolute error CDF", xlab="",ylab="")
		mtext("correct to within (hrs)",side=1,line=2.2,cex=0.8)
		mtext(paste("% correct (N = ",length(hrerr),")",sep=""),side=2,line=2.4,cex=0.8)
		abline(h=c(50,80,100),v=c(median(abs(hrerr)),quantile(abs(hrerr),0.8)),col="grey")
		asTime <- function(x){
			hrs <- (x%%24)%/%1; 
			min <- round(60*(x%%1)); 
			if(min==60){hrs <- (hrs+1)%%24; min<-0}
			sprintf("%2i:%02i",hrs,min)
		}
		text(quantile(abs(hrerr),0.8),30,asTime(quantile(abs(hrerr),0.8)),cex=0.9)
		text(median(abs(hrerr)),10, asTime(median(abs(hrerr))),cex=0.9)
	}
	points(hrsoff, 100-fracacc, col=col)
	lines(hrsoff, 100-fracacc, col=col)
	norm.fracacc <- (100-fracacc)/100
	norm.hrsoff <- hrsoff/12
	auc <- sum(norm.fracacc[-1]*diff(norm.hrsoff))
	if(!add){
		text(10,20, sprintf("nAUC=%.2f",auc),cex=1,font=2)
	}
	invisible(auc)
}


predplot <- function(trueHr,predHr,col=1,pch=1,main="Time of Day (24h)",...){
	# do both of the above -- set par(mfrow=c(2,1)) or (1,2) first!
	opar <- par(xpd=F,mar=c(4,4,3,1))	
	on.exit(par(opar))
	out <- timeErr(trueHr,predHr,plot=T,col,pch,main) 
	tolplot(trueHr,predHr)
	invisible(out)
}