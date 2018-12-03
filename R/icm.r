#icm
icm.ft=function(x, d, wt){

  #1. evt: # of death at x_m; tot: # of sub inspected at x_m 
  x.unq=unique(x)
  M=length(x.unq);
  
  evt=tot=rep(NA,M) 
  for(m in 1:M){
    evt[m]=sum((x.unq[m]==x & d==1)*wt)
    tot[m]=sum((x.unq[m]==x)*wt)
  }
  
  #2. remove zero weight & pava
  idx.zwt=which(tot>0)
  
  if(length(idx.zwt)==0){
    F.zwt=x.unq.zwt=NA
  }else{
    evt.zwt=evt[idx.zwt]
    tot.zwt=tot[idx.zwt]
    x.unq.zwt=x.unq[idx.zwt]

    F.zwt=pava((evt.zwt/tot.zwt), tot.zwt)
  }

  return(list(F=F.zwt, x=x.unq.zwt))
}
