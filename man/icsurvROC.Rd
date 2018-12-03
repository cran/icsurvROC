\name{icsurvROC}
\alias{icsurvROC}
\title{Interval Censored survival ROC}
\description{Nonparametric and semiparametric estimations of the time-dependent ROC curve for the interval-censored failure time data.}
\usage{
  icsurvROC(Time, Status, Marker, pred.time, method, wt=NULL, span=NULL)
}
\arguments{
  \item{Time}{Monitoring time}
  \item{Status}{Event indicator (1: event occurs before or at monitoring time; 0: event occurs after monitoring time)}
  \item{Marker}{Predictior or marekr value}
  \item{pred.time}{Prediction time of the ROC curve}
  \item{method}{"np" for nonparametric model using local NPMLE; "sp" for semiparametric proportional hazard model}
  \item{wt}{Weight, such as inverse probablity weighting}
  \item{span}{Smoothing bandwidth parameter for method="np"}
}
\details{It estimates a time-dependent ROC curve for the current status data based on the local NPMLE if method="np" or proportional hazards model if method="sp". For method="np", optimal bandwidth selection based on maximum likelihood cross varlidation is used if span is not specified or span=NUL. For method="sp", span is not needed to be specified.}

\references{Yunro Chung, Tianxi Cai, Yingye Zheng, Estimating Diagnostic Accuracy Measures for Current Status Survival Data with Application to Prostate Cancer Active Surveillance Study (in progress)}
\author{Yunro Chung [cre]}

\examples{
Time=   c(1,2,5,3,9,8,9,4,6,4)
Status= c(1,1,1,0,1,1,1,0,0,0)
Marker= c(8,2,6,3,1,4,5,1,3,7)

#np at year 3
nobs=length(Time)
span=sd(Marker)*nobs^(-1/7)
RES1=icsurvROC(Time, Status, Marker, pred.time=3, method="np", span=span)
print(RES1)

#sp at year 3
RES2=icsurvROC(Time, Status, Marker, pred.time=3, method="sp")
print(RES2)
}

\keyword{Semiparametric Cox proportional hazards, local NPMLE, maximum likelihood croos validation}
