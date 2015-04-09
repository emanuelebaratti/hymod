

################################################################
### Hymod: R-language
################################################################

hymodr=function(param,area,tdelta,e,p,w_initial=0,wslow_initial=0,wquick_initial=0)
{
  
  # parameters #
  cmax=param[1]
  beta=param[2]
  alfa=param[3]
  kslow=param[4]
  kquick=param[5]
  
  # Wathershed area: conversion from km^2 to m^2 
  area=area*1000000
  
  # Set the conversion factor
  fatconv=1/1000/tdelta*area
  
  # Initialisation of variables
  ntstep=length(p)
  w2<-w_initial
  wslow<-wslow_initial
  wquick<-rep(wquick_initial,3) 
  c1<-0
  w1<-rep(0,ntstep)
  er1<-rep(0,ntstep)
  er2<-rep(0,ntstep)
  er<-rep(0,ntstep)
  ep<-rep(0,ntstep)
  qslow<-rep(0,ntstep)
  qtquick<-rep(0,ntstep)
  qtslow<-rep(0,ntstep)
  qtot<-rep(0,ntstep)
  # Computation loop 
  for ( i in 1:ntstep) {
    w1[i]<-w2
      dummy<-(1-((beta+1)*w1[i]/cmax))
      dummy<-max(dummy,0)
    c1<-cmax*(1-(dummy^(1/(beta +1))))
      rm(dummy)
    c2<-min(c1+p[i],cmax)
    er1[i]<-max((p[i]-cmax+c1),0)
      dummy<-1-c2/cmax
      dummy<-max(dummy,0)
    w2<-(cmax/(beta+1))*(1-(dummy^(beta+1)))  
      rm(dummy)
    er2[i]<-max((c2-c1)-(w2-w1[i]),0)
    ep[i]<-(1-(((cmax-c2)/(beta+1))/(cmax/(beta+1))))*e[i]
    w2<-max(w2-ep[i],0)
    
    # Subdivision of the surface runoff
    uquick<-alfa*er2[i]+er1[i]
    uslow<-(1-alfa)*er2[i]
    
    # Slow flow component
    wslow<-(1-kslow)*wslow+(1-kslow)*uslow
    qslow[i]<-(kslow/(1-kslow))*wslow
    
    # Quick flow component
    qquick<-0 
    for ( j in 1:3)
    {
      wquick[j]<-(1-kquick)*wquick[j]+(1-kquick)*uquick
      qquick<-(kquick/(1-kquick))*wquick[j]
      uquick<-qquick
    }
    # Quick, Slow and Total surface runoff
    qtslow[i]<-qslow[i]*fatconv
    qtquick[i]<-qquick*fatconv
    qtot[i]<-qtquick[i]+qtslow[i]  
    
  }
  
  return(list(w_tot=w1,ep_tot=ep,q_slow=qtslow,q_quick=qtquick,q_tot=qtot))
  
}


################################################################
### Hymod: FORTRAN 77
################################################################

hymodf<-function(param,
                 area,
                 tdelta,
                 e,
                 p,
                 w_initial=0,
                 wslow_initial=0,
                 wquick_initial=0)  {
  nstep<-length(p)
  qout<-rep(999,nstep)
  qsimf<-.Fortran("hymodfortran", 
                    param=as.double(param),
                    area=as.double(area),
                    tdelta=as.integer(tdelta),
                    e=as.double(e),
                    p=as.double(p),
                    nstep=as.integer(nstep),
                    winitial=as.double(w_initial),
                    wslowinitial=as.double(wslow_initial),
                    wquickinitial=as.double(wquick_initial),
                    w1=as.double(qout),
                    evapt=as.double(qout),
                    qtslow=as.double(qout),
                    qtquick=as.double(qout),
                    qtot=as.double(qout),
                    PACKAGE="hymod")
  return(list(w_tot=qsimf$w1,ep_tot=qsimf$evapt,
              q_slow=qsimf$qtslow,q_quick=qsimf$qtquick,q_tot=qsimf$qtot))
  
}
