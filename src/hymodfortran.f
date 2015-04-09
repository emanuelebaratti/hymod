      subroutine hymodfortran(param,area,tdelta,e,p,
     *     ntstep,winitial,wslowinitial,
     *	   wquickinitial,w1,evapt,
     *     qtslow,qtquick,qtot)
!
!  Declarations for subroutine arguments.
!
      implicit none
      integer tdelta,ntstep,deltat,i,j
      double precision er1,er2,
     *     e(ntstep),p(ntstep),
     *     winitial,wslowinitial,wquickinitial,
     *     qtquick(ntstep),qtslow(ntstep),
     *     qtot(ntstep),evapt(ntstep),param(5),
     *	   area,cmax,betap,alfa,kslow,kquick,
     *     w1(ntstep),w2,fatconv,c1,c2,qinitial,
     *     wquick(3),qquick,wslow,uslow,qslow,
     *     uquick,dummy


!  Parameters
	cmax=param(1)
	betap=param(2)
	alfa=param(3)
	kslow=param(4)
	kquick=param(5)

!  Wathershed area: conversion from kmq to mq
	area=area*1000000.

!  Set the conversion factor
	fatconv=(1./1000.0/tdelta)*area 
  
!  Store initialization
	w1=0.0 
	w2=winitial 
	c1=0.0 
	c2=0.0 
	uquick=0.0 
	uslow=0.0 
 	wslow=qinitial/(kslow*fatconv)
	 do i=1,3
    	   wquick(i)=0.0
         end do

!---------------------------------------------------
!  Model loop  
!
	 do i=1,ntstep
!  Compute excess precipitation and evaporation
    	   w1(i)=w2
           ! remove negative number
           	dummy=(1. - ((betap + 1.) * w1(i)/cmax))
        	dummy=max(dummy,0.0)
           c1=cmax * (1. - (dummy**(1./(betap + 1.))))
	   c2=min(c1+p(i),cmax)  
	   er1=max((p(i) - cmax + c1),0.0)
	   	dummy=1. - c2/cmax
	   	dummy=max(dummy,0.0)
	   w2=(cmax/(betap + 1.))*(1. - (dummy**(betap + 1.)))
	   er2=max((c2 - c1) - (w2 - w1(i)),0.0)
	   evapt(i)=(1.-(((cmax-c2)/(betap+ 1.))/(cmax/(betap+ 1.))))*e(i)
	   w2=max(w2 - evapt(i),0.0)  
    
!  Partition uquick and uslow into quick and slow flow component 
	   uquick=alfa*er2+er1
	   uslow=(1.- alfa)*er2

!  Route slow flow component with single linear reservoir (kslow)
	   wslow=(1. - kslow)*wslow + (1.- kslow)*uslow
	   qslow=(kslow/(1. - kslow))*wslow
      
 	   qquick=0.0
	     do j=1,3
		wquick(j)=(1.- kquick)*wquick(j) + (1. - kquick)*uquick
		qquick=(kquick/(1.- kquick))*wquick(j)
		uquick=qquick
 	     end do
	    
!  Compute quick, slow and total flow for timestep
 	   qtslow(i)=qslow*fatconv
	   qtquick(i)=qquick*fatconv
  	   qtot(i)=qtquick(i) + qtslow(i)
  
	 end do

!---------------------------------------------------

      end subroutine hymodfortran

