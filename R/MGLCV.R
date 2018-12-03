MGLVC.ft=function(data.k1.ord,res.cv,span,n){
  x=data.k1.ord$X
  d=data.k1.ord$D
  y=data.k1.ord$Y
  ipw=data.k1.ord$ipw
  
  y.cv=res.cv$y
  
  x.unq=unique(x)
  M=length(x.unq);
  nk=length(x) #sample for kth fold
  lik=rep(NA,nk)
  for(i in 1:nk){
    #(1) compute wt evt and tot
    y.span=(y-y[i])/span
    num.wt=dnorm(y.span)
    den.wt=sum(num.wt)
    wt=num.wt/den.wt
    dwt=wt*ipw #double weight
    
    evt=tot=rep(NA,M) 
    for(m in 1:M){
      evt[m]=sum((x.unq[m]==x & d==1)*dwt)
      tot[m]=sum((x.unq[m]==x)*dwt)
    }
    
    #(2) nearest point: j(i)
    ji=which.min(abs(y[i]-y.cv))
    res.ji=res.cv[[ji]]
    
    F.ji=res.ji$F
    x.ji=res.ji$x

    #(3) F.ji is defined on x.unq by right-cont assumption
    F.tilder=rep(0,M)
    if(x.unq[1]<x.ji[1]){
      m.idx1=max(which(x.unq<x.ji[1]))
      F.tilder[1:m.idx1]=0
      for(m in (m.idx1+1):M)
        F.tilder[m]=F.ji[which.max(which(x.ji<=x.unq[m]))]
    }else{
      for(m in 1:M)
        F.tilder[m]=F.ji[which.max(which(x.ji<=x.unq[m]))]
    }    

    #(4) compute log lik
    idx=which(0<F.tilder & F.tilder<1) #to prevent lik is infinite
    lik[i]=sum((evt*log(F.tilder)+(tot-evt)*log(1-F.tilder))[idx]) #lik
  }
  lik2=sum(lik)
  
  return(lik2)
}