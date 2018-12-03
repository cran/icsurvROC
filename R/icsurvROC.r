icsurvROC=function(Time, Status, Marker, pred.time, method, wt=NULL, span=NULL){
  res=icROC(X=Time,D=Status,Y=Marker, wt=wt, eval.time=pred.time, span=span, method=method)
  plot(res$tpf~res$fpf,ylab='TPF',xlab='FPF',main=paste("ROC at t=",pred.time,sep=""),type='l',lwd=1.5)
  abline(0,1,col='gray', lwd=1.5)
  return(res)
}
