icROC=function(X,D,Y, wt=NULL, eval.time, span=NULL, method){
  tpf=fpf=auc=y=NA

  if(sum(is.na(X+D+Y+wt))){
    print("NA is not allowed")
  }else{
    ###
    #1. sorting
    ###
    n=length(X)
    if(is.null(wt))
      wt=rep(1,n)

    data=data.frame(X=X,D=D,Y=Y,ipw=wt)
    data.ord=data[order(data$X),]

    x=data.ord$X
    d=data.ord$D
    y=data.ord$Y
    ipw=data.ord$ipw
    ####
    #2. Compute St.y
    ###
    St.y=rep(NA,n)
    if(method=="sp"){ #IC cox
      data.ic=cs2ic(x,as.logical(d))
      data.ic.Surv=Surv(time=data.ic[,1], time2=data.ic[,2], type='interval2')
      options(warn=2)
      cox.warning=try(ph.fit<-ic_sp(data.ic.Surv ~ y, model="ph",weights=ipw),silent=TRUE)
      options(warn=0)

      if(class(cox.warning)=="try-error"){
        print("Proportional hazards model does not converge")
      }else{
        ph.fit.S=getSCurves(ph.fit)

        bS=ph.fit.S$S_curves$baseline #baseline hazard
        bt=ph.fit.S$Tbull_ints[, 2]   #time, x-axis
        beta=unname(ph.fit$coefficients)

        eval.pt=tail(which(bt<=eval.time),1)
        St.y=bS[eval.pt]^(exp(y*beta))
      }
    }else if(method=="np"){ #Local NPMLE
      if(is.null(span)){ #K=5-fold ML cross validation for choosing h_opt
        k.fold=5
        cv.ind=runif(n,0,1)
        cv.idx=list()
        cv.idx[[1]]=which(cv.ind<=0.2)
        cv.idx[[2]]=which(0.2<cv.ind & cv.ind<=0.4)
        cv.idx[[3]]=which(0.4<cv.ind & cv.ind<=0.6)
        cv.idx[[4]]=which(0.6<cv.ind & cv.ind<=0.8)
        cv.idx[[5]]=which(0.8<cv.ind)

        span.all=seq(0.1,3,0.1)*sd(y)  #candidate
        n.span=length(span.all)
        GMLCV=rep(NA,n.span)
        for(j in 1:n.span){
          lik=rep(NA,k.fold)
          for(k in 1:k.fold){
            cv.idx1=cv.idx[[k]] #(a) D^k: k fold data
            cv.idx2=cv.idx      #(b) D^-k: data excluding k fold
            cv.idx2[[k]]=NULL
            cv.idx2=unlist(cv.idx2)

            data.k1=data[cv.idx1,] #(a) D^k
            data.k2=data[cv.idx2,] #(b) D^-k

            data.k1.ord=data.k1[order(data.k1$X),] #(a) D^k, ordered by X
            data.k2.ord=data.k2[order(data.k2$X),] #(b) D^-k, ordered by X

            res.cv=list()
            for(i in 1:nrow(data.k2)){ #using D^-k
              #(1) NW weights
              y.span=(data.k2.ord$Y-data.k2.ord$Y[i])/span.all[j]
              num.wt=dnorm(y.span)
              den.wt=sum(num.wt)
              wt=num.wt/den.wt
              dwt=wt*data.k2.ord$ipw #double weight

              #(2) localized NPMLE
              res.cv[[i]]=icm.ft(data.k2.ord$X,data.k2.ord$D,dwt)
            }
            res.cv$y=data.k2.ord$Y
            lik[k]=MGLVC.ft(data.k1.ord,res.cv,span.all[j],n) #lik for kth folder with jth span
          }
          GMLCV[j]=sum(lik)
        }
        span=span.all[which.max(GMLCV)]
        #plot(GMLCV~span.all)
        #print(span)
      }

      for(i in 1:n){ #Local NPMLE
        #(1) NW weights
        y.span=(y-y[i])/span
        num.wt=dnorm(y.span)
        den.wt=sum(num.wt)
        wt=num.wt/den.wt
        dwt=wt*ipw #double weight

        #(2) localized NPMLE
        res=icm.ft(x,d,dwt)
        eval.pt=tail(which(res$x<=eval.time),1) #right cont
        if(length(eval.pt)==0) eval.pt=1
        F.y=res$F[eval.pt][1]
        St.y[i]=1-F.y
      }
    }

    ###
    #4. mzr
    ###
    #4.1.initial
    n.ipw=sum(ipw)
    y.unq=sort(unique(y))
    St=sum(ipw*St.y)/n.ipw

    #4.2 tpf&fpf at all y.unq
    n.unq=length(y.unq)
    tpf.a=fpf.a=rep(NA,n.unq)
    for(i in 1:n.unq){
      Sct=sum(ipw*St.y*(y>y.unq[i]))/n.ipw
      Fc=sum(ipw*(y<=y.unq[i]))/n.ipw
      if(St!=1) tpf.a[i]=(1-Fc-Sct)/(1-St)
      if(St!=0) fpf.a[i]=Sct/St
    }
    tpf.a[tpf.a>1]=1; tpf.a[tpf.a<0]=0;
    fpf.a[fpf.a>1]=1; fpf.a[fpf.a<0]=0;

    #4.3 sorted by y
    data0.a=data.frame(tpf.a,fpf.a,y.unq)
    data0.a=data0.a[order(data0.a$y),]
    y.a=data0.a$y
    tpf.a=data0.a$tpf.a
    fpf.a=data0.a$fpf.a

    #4.4 including the end points
    y.a=c(-Inf,y.a,Inf)
    tpf.a=c(1,tpf.a,0)
    fpf.a=c(1,fpf.a,0)
    n.a=n.unq+2

    #4.5 mid pt integration, like the same approach as SurvivalROC package, Heagerty et al.
    d.fpf=fpf.a[-n.a]-fpf.a[-1]
    mid.tpf=(tpf.a[-n.a] + tpf.a[-1])/2
    auc=sum(d.fpf*mid.tpf)
  }
  return(list(tpf=tpf.a,
              fpf=fpf.a,
              cut.value=y.a,
              auc=auc,
              eval.time=eval.time))
}

