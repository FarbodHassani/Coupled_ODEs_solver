      subroutine dopri5(n,fcn,x,y,xend,                                 &
     &     rtol,atol,itol,                                              &
     &     solout,iout,                                                 &
     &     work,lwork,iwork,liwork,idid)                      
! ----------------------------------------------------------            
!     numerical solution of a system of first 0rder                     
!     ordinary differential equations  y'=f(x,y).                       
!     this is an explicit runge-kutta method of order (4)5              
!     due to dormand & prince (with stepsize control and                
!     dense output).                                                    
!                                                                       
!     authors: e. hairer and g. wanner                                  
!              universite de geneve, dept. de mathematiques             
!              ch-1211 geneve 24, switzerland                           
!              e-mail:  ernst.hairer@math.unige.ch                      
!                       gerhard.wanner@math.unige.ch                    
!                                                                       
!     this code is described in:                                        
!         e. hairer, s.p. norsett and g. wanner, solving ordinary       
!         differential equations i. nonstiff problems. 2nd edition.     
!         springer series in computational mathematics,                 
!         springer-verlag (1993)                                        
!                                                                       
!     version of april 25, 1996                                         
!     (latest correction of a small bug: august 8, 2005)                
!                                                                       
!     input parameters                                                  
!     ----------------                                                  
!     n           dimension of the system                               
!                                                                       
!     fcn         name (external) of subroutine computing the           
!                 value of f(x,y):                                      
!                    subroutine fcn(n,x,y,f)                  
!                    double precision x,y(n),f(n)                       
!                    f(1)=...   etc.                                    
!                                                                       
!     x           initial x-value                                       
!                                                                       
!     y(n)        initial values for y                                  
!                                                                       
!     xend        final x-value (xend-x may be positive or negative)    
!                                                                       
!     rtol,atol   relative and absolute error tolerances. they          
!                 can be both scalars or else both vectors of length n. 
!                                                                       
!     itol        switch for rtol and atol:                             
!                   itol=0: both rtol and atol are scalars.             
!                     the code keeps, roughly, the local error of       
!                     y(i) below rtol*abs(y(i))+atol                    
!                   itol=1: both rtol and atol are vectors.             
!                     the code keeps the local error of y(i) below      
!                     rtol(i)*abs(y(i))+atol(i).                        
!                                                                       
!     solout      name (external) of subroutine providing the           
!                 numerical solution during integration.                
!                 if iout.ge.1, it is called after every successful step
!                 supply a dummy subroutine if iout=0.                  
!                 it must have the form                                 
!                    subroutine solout (nr,xold,x,y,n,con,icomp,nd,     
!                                       ipar,irtrn)                
!                    dimension y(n),con(5*nd),icomp(nd)                 
!                    ....                                               
!                 solout furnishes the solution "y" at the nr-th        
!                    grid-point "x" (thereby the initial value is       
!                    the first grid-point).                             
!                 "xold" is the preceeding grid-point.                  
!                 "irtrn" serves to interrupt the integration. if irtrn 
!                    is set <0, dopri5 will return to the calling progra
!                    if the numerical solution is altered in solout,    
!                    set  irtrn = 2                                     
!                                                                       
!          -----  continuous output: -----                              
!                 during calls to "solout", a continuous solution       
!                 for the interval [xold,x] is available through        
!                 the function                                          
!                        >>>   contd5(i,s,con,icomp,nd)   <<<           
!                 which provides an approximation to the i-th           
!                 component of the solution at the point s. the value   
!                 s should lie in the interval [xold,x].                
!                                                                       
!     iout        switch for calling the subroutine solout:             
!                    iout=0: subroutine is never called                 
!                    iout=1: subroutine is used for output.             
!                    iout=2: dense output is performed in solout        
!                            (in this case work(5) must be specified)   
!                                                                       
!     work        array of working space of length "lwork".             
!                 work(1),...,work(20) serve as parameters for the code.
!                 for standard use, set them to zero before calling.    
!                 "lwork" must be at least  8*n+5*nrdens+21             
!                 where  nrdens = iwork(5)                              
!                                                                       
!     lwork       declared lenght of array "work".                      
!                                                                       
!     iwork       integer working space of lenght "liwork".             
!                 iwork(1),...,iwork(20) serve as parameters for the cod
!                 for standard use, set them to zero before calling.    
!                 "liwork" must be at least nrdens+21 .                 
!                                                                       
!     liwork      declared lenght of array "iwork".                     
!  
!-----------------------------------------------------------------------
!                                                                       
!     sophisticated setting of parameters                               
!     -----------------------------------                               
!              several parameters (work(1),...,iwork(1),...) allow      
!              to adapt the code to the problem and to the needs of     
!              the user. for zero input, the code chooses default values
!                                                                       
!    work(1)   uround, the rounding unit, default 2.3d-16.              
!                                                                       
!    work(2)   the safety factor in step size prediction,               
!              default 0.9d0.                                           
!                                                                       
!    work(3), work(4)   parameters for step size selection              
!              the new step size is chosen subject to the restriction   
!                 work(3) <= hnew/hold <= work(4)                       
!              default values: work(3)=0.2d0, work(4)=10.d0             
!                                                                       
!    work(5)   is the "beta" for stabilized step size control           
!              (see section iv.2). larger values of beta ( <= 0.1 )     
!              make the step size control more stable. dopri5 needs     
!              a larger beta than higham & hall. negative work(5)       
!              provoke beta=0.                                          
!              default 0.04d0.                                          
!                                                                       
!    work(6)   maximal step size, default xend-x.                       
!                                                                       
!    work(7)   initial step size, for work(7)=0.d0 an initial guess     
!              is computed with help of the function hinit              
!                                                                       
!    iwork(1)  this is the maximal number of allowed steps.             
!              the default value (for iwork(1)=0) is 100000.            
!                                                                       
!    iwork(2)  switch for the choice of the coefficients                
!              if iwork(2).eq.1  method dopri5 of dormand and prince    
!              (table 5.2 of section ii.5).                             
!              at the moment this is the only possible choice.          
!              the default value (for iwork(2)=0) is iwork(2)=1.        
!                                                                       
!    iwork(3)  switch for printing error messages                       
!              if iwork(3).lt.0 no messages are being printed           
!              if iwork(3).gt.0 messages are printed with               
!              write (iwork(3),*) ...                                   
!              default value (for iwork(3)=0) is iwork(3)=6             
!                                                                       
!    iwork(4)  test for stiffness is activated after step number        
!              j*iwork(4) (j integer), provided iwork(4).gt.0.          
!              for negative iwork(4) the stiffness test is              
!              never activated; default value is iwork(4)=1000          
!                                                                       
!    iwork(5)  = nrdens = number of components, for which dense output  
!              is required; default value is iwork(5)=0;                
!              for   0 < nrdens < n   the components (for which dense   
!              output is required) have to be specified in              
!              iwork(21),...,iwork(nrdens+20);                          
!              for  nrdens=n  this is done by the code.                 
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     output parameters                                                 
!     -----------------                                                 
!     x           x-value for which the solution has been computed      
!                 (after successful return x=xend).                     
!                                                                       
!     y(n)        numerical solution at x                               
!                                                                       
!     h           predicted step size of the last accepted step         
!                                                                       
!     idid        reports on successfulness upon return:                
!                   idid= 1  computation successful,                    
!                   idid= 2  comput. successful (interrupted by solout) 
!                   idid=-1  input is not consistent,                   
!                   idid=-2  larger nmax is needed,                     
!                   idid=-3  step size becomes too small.               
!                   idid=-4  problem is probably stiff (interrupted).   
!                                                                       
!   iwork(17)  nfcn    number of function evaluations                   
!   iwork(18)  nstep   number of computed steps                         
!   iwork(19)  naccpt  number of accepted steps                         
!   iwork(20)  nrejct  number of rejected steps (due to error test),    
!                      (step rejections in the first step are not counte
!-----------------------------------------------------------------------
! *** *** *** *** *** *** *** *** *** *** *** *** ***                   
!          declarations                                                 
! *** *** *** *** *** *** *** *** *** *** *** *** ***                   
      implicit double precision (a-h,o-z) 
      dimension y(n),atol(*),rtol(*),work(lwork),iwork(liwork) 
      logical arret 
      external fcn,solout 
! *** *** *** *** *** *** ***                                           
!        setting the parameters                                         
! *** *** *** *** *** *** ***                                           
      nfcn=0 
      nstep=0 
      naccpt=0 
      nrejct=0 
      arret=.false. 
! -------- iprint for monitoring the printing                           
      if(iwork(3).eq.0)then 
         iprint=6 
      else 
         iprint=iwork(3) 
      end if 
! -------- nmax , the maximal number of steps -----                     
      if(iwork(1).eq.0)then 
         nmax=100000 
      else 
         nmax=iwork(1) 
         if(nmax.le.0)then 
            if (iprint.gt.0) write(iprint,*)                            &
     &           ' wrong input iwork(1)=',iwork(1)                      
            arret=.true. 
         end if 
      end if 
! -------- meth   coefficients of the method                            
      if(iwork(2).eq.0)then 
         meth=1 
      else 
         meth=iwork(2) 
         if(meth.le.0.or.meth.ge.4)then 
            if (iprint.gt.0) write(iprint,*)                            &
     &           ' curious input iwork(2)=',iwork(2)                    
            arret=.true. 
         end if 
      end if 
! -------- nstiff   parameter for stiffness detection                   
      nstiff=iwork(4) 
      if (nstiff.eq.0) nstiff=1000 
      if (nstiff.lt.0) nstiff=nmax+10 
! -------- nrdens   number of dense output components                   
      nrdens=iwork(5) 
      if(nrdens.lt.0.or.nrdens.gt.n)then 
         if (iprint.gt.0) write(iprint,*)                               &
     &        ' curious input iwork(5)=',iwork(5)                       
         arret=.true. 
      else 
         if(nrdens.gt.0.and.iout.lt.2)then 
            if (iprint.gt.0) write(iprint,*)                            &
     &           ' warning: put iout=2 for dense output '               
         end if 
         if (nrdens.eq.n) then 
            do  i=1,nrdens 
               iwork(20+i)=i 
            end do 
         end if 
      end if 
! -------- uround   smallest number satisfying 1.d0+uround>1.d0         
      if(work(1).eq.0.d0)then 
         uround=2.3d-16 
      else 
         uround=work(1) 
         if(uround.le.1.d-35.or.uround.ge.1.d0)then 
            if (iprint.gt.0) write(iprint,*)                            &
     &           ' which machine do you have? your uround was:',work(1) 
            arret=.true. 
         end if 
      end if 
! -------  safety factor -------------                                  
      if(work(2).eq.0.d0)then 
         safe=0.9d0 
      else 
         safe=work(2) 
         if(safe.ge.1.d0.or.safe.le.1.d-4)then 
            if (iprint.gt.0) write(iprint,*)                            &
     &           ' curious input for safety factor work(2)=',work(2)    
            arret=.true. 
         end if 
      end if 
! -------  fac1,fac2     parameters for step size selection             
      if(work(3).eq.0.d0)then 
         fac1=0.2d0 
      else 
         fac1=work(3) 
      end if 
      if(work(4).eq.0.d0)then 
         fac2=10.d0 
      else 
         fac2=work(4) 
      end if 
! --------- beta for step control stabilization -----------             
      if(work(5).eq.0.d0)then 
         beta=0.04d0 
      else 
         if(work(5).lt.0.d0)then 
            beta=0.d0 
         else 
            beta=work(5) 
            if(beta.gt.0.2d0)then 
               if (iprint.gt.0) write(iprint,*)                         &
     &              ' curious input for beta: work(5)=',work(5)         
               arret=.true. 
            end if 
         end if 
      end if 
! -------- maximal step size                                            
      if(work(6).eq.0.d0)then 
         hmax=xend-x 
      else 
         hmax=work(6) 
      end if 
! -------- initial step size                                            
      h=work(7) 
! ------- prepare the entry-points for the arrays in work -----         
      iey1=21 
      iek1=iey1+n 
      iek2=iek1+n 
      iek3=iek2+n 
      iek4=iek3+n 
      iek5=iek4+n 
      iek6=iek5+n 
      ieys=iek6+n 
      ieco=ieys+n 
! ------ total storage requirement -----------                          
      istore=ieys+5*nrdens-1 
      if(istore.gt.lwork)then 
         if (iprint.gt.0) write(iprint,*)                               &
     &        ' insufficient storage for work, min. lwork=',istore      
         arret=.true. 
      end if 
      icomp=21 
      istore=icomp+nrdens-1 
      if(istore.gt.liwork)then 
         if (iprint.gt.0) write(iprint,*)                               &
     &        ' insufficient storage for iwork, min. liwork=',istore    
         arret=.true. 
      end if 
! ------ when a fail has occured, we return with idid=-1                
      if (arret) then 
         idid=-1 
         return 
      end if 
! -------- call to core integrator ------------                         
      call dopcor(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,iprint,          &
     &     solout,iout,idid,nmax,uround,meth,nstiff,safe,beta,fac1,fac2,&
     &     work(iey1),work(iek1),work(iek2),work(iek3),work(iek4),      &
     &     work(iek5),work(iek6),work(ieys),work(ieco),iwork(icomp),    &
     &     nrdens,nfcn,nstep,naccpt,nrejct)                   
      work(7)=h 
      iwork(17)=nfcn 
      iwork(18)=nstep 
      iwork(19)=naccpt 
      iwork(20)=nrejct 
! ----------- return -----------                                        
      return 
      END                                           
!                                                                       
!                                                                       
!                                                                       
!  ----- ... and here is the core integrator  ----------                
!                                                                       
      subroutine dopcor(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,iprint,    &
     &     solout,iout,idid,nmax,uround,meth,nstiff,safe,beta,fac1,fac2,&
     &     y1,k1,k2,k3,k4,k5,k6,ysti,cont,icomp,nrd,          &
     &     nfcn,nstep,naccpt,nrejct)                                    
! ----------------------------------------------------------            
!     core integrator for dopri5                                        
!     parameters same as in dopri5 with workspace added                 
! ----------------------------------------------------------            
!         declarations                                                  
! ----------------------------------------------------------            
      implicit double precision (a-h,o-z) 
      double precision k1(n),k2(n),k3(n),k4(n),k5(n),k6(n) 
      dimension y(n),y1(n),ysti(n),atol(*),rtol(*) 
      dimension cont(5*nrd),icomp(nrd) 
      logical reject,last 
      external fcn 
      common /condo5/xold,hout
! *** *** *** *** *** *** ***                                           
!  initialisations                                                      
! *** *** *** *** *** *** ***                                           
      if (meth.eq.1) call cdopri(c2,c3,c4,c5,e1,e3,e4,e5,e6,e7,         &
     &     a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,                     &
     &     a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,                     &
     &     d1,d3,d4,d5,d6,d7)                                           
      facold=1.d-4 
      expo1=0.2d0-beta*0.75d0 
      facc1=1.d0/fac1 
      facc2=1.d0/fac2 
      posneg=sign(1.d0,xend-x) 
! --- initial preparations                                              
      atoli=atol(1) 
      rtoli=rtol(1) 
      last=.false. 
      hlamb=0.d0 
      iasti=0 
      call fcn(n,x,y,k1) 
      hmax=abs(hmax) 
      iord=5 
      if (h.eq.0.d0) h=hinit(n,fcn,x,y,xend,posneg,k1,k2,k3,iord,       &
     &     hmax,atol,rtol,itol)                               
      nfcn=nfcn+2 
      reject=.false. 
      xold=x 
      if (iout.ne.0) then 
         irtrn=1 
         hout=h

         call solout(naccpt+1,xold,x,y,n,cont,icomp,nrd,                &
              &        irtrn)
         if (irtrn.lt.0) goto 79 
      else 
         irtrn=0 
      end if 
! --- basic integration step                                            
    1 continue 
      if (nstep.gt.nmax) goto 78 
      if (0.1d0*abs(h).le.abs(x)*uround)goto 77 
      if ((x+1.01d0*h-xend)*posneg.gt.0.d0) then 
         h=xend-x 
         last=.true. 
      end if 
      nstep=nstep+1 
! --- the first 6 stages                                                
      if (irtrn.ge.2) then 
         call fcn(n,x,y,k1) 
      end if 
      do  i=1,n 
         y1(i)=y(i)+h*a21*k1(i) 
      end do 
      call fcn(n,x+c2*h,y1,k2) 
      do  i=1,n 
         y1(i)=y(i)+h*(a31*k1(i)+a32*k2(i)) 
      end do 
      call fcn(n,x+c3*h,y1,k3) 
      do i=1,n 
         y1(i)=y(i)+h*(a41*k1(i)+a42*k2(i)+a43*k3(i)) 
      end do 
      call fcn(n,x+c4*h,y1,k4) 
      do i=1,n 
         y1(i)=y(i)+h*(a51*k1(i)+a52*k2(i)+a53*k3(i)+a54*k4(i)) 
      end do 
      call fcn(n,x+c5*h,y1,k5) 
      do  i=1,n 
      ysti(i)=y(i)+h*(a61*k1(i)+a62*k2(i)+a63*k3(i)+a64*k4(i)+a65*k5(i)) 
      end do 
      xph=x+h 
      call fcn(n,xph,ysti,k6) 
      do  i=1,n 
      y1(i)=y(i)+h*(a71*k1(i)+a73*k3(i)+a74*k4(i)+a75*k5(i)+a76*k6(i)) 
      end do 
      call fcn(n,xph,y1,k2) 
      if (iout.ge.2) then 
         do 40 j=1,nrd 
            i=icomp(j) 
            cont(4*nrd+j)=h*(d1*k1(i)+d3*k3(i)+d4*k4(i)+d5*k5(i)        &
     &           +d6*k6(i)+d7*k2(i))                                    
   40    continue 
      end if 
      do  i=1,n 
         k4(i)=(e1*k1(i)+e3*k3(i)+e4*k4(i)+e5*k5(i)+e6*k6(i)+e7*k2(i))*h 
      end do 
      nfcn=nfcn+6 
! --- error estimation                                                  
      err=0.d0 
      if (itol.eq.0) then 
         do  i=1,n 
            sk=atoli+rtoli*max(abs(y(i)),abs(y1(i))) 
            err=err+(k4(i)/sk)**2 
         end do 
      else 
         do  i=1,n 
            sk=atol(i)+rtol(i)*max(abs(y(i)),abs(y1(i))) 
            err=err+(k4(i)/sk)**2 
         end do 
      end if 
      err=sqrt(err/n) 
! --- computation of hnew                                               
      fac11=err**expo1 
! --- lund-stabilization                                                
      fac=fac11/facold**beta 
! --- we require  fac1 <= hnew/h <= fac2                                
      fac=max(facc2,min(facc1,fac/safe)) 
      hnew=h/fac 
      if(err.le.1.d0)then 
! --- step is accepted                                                  
         facold=max(err,1.0d-4) 
         naccpt=naccpt+1 
! ------- stiffness detection                                           
         if (mod(naccpt,nstiff).eq.0.or.iasti.gt.0) then 
            stnum=0.d0 
            stden=0.d0 
            do 64 i=1,n 
               stnum=stnum+(k2(i)-k6(i))**2 
               stden=stden+(y1(i)-ysti(i))**2 
   64       continue 
            if (stden.gt.0.d0) hlamb=h*sqrt(stnum/stden) 
            if (hlamb.gt.3.25d0) then 
               nonsti=0 
               iasti=iasti+1 
               if (iasti.eq.15) then 
                  if (iprint.gt.0) write (iprint,*)                     &
     &                 ' the problem seems to become stiff at x = ',x   
                  if (iprint.le.0) goto 76 
               end if 
            else 
               nonsti=nonsti+1 
               if (nonsti.eq.6) iasti=0 
            end if 
         end if 
         if (iout.ge.2) then 
            do 43 j=1,nrd 
               i=icomp(j) 
               yd0=y(i) 
               ydiff=y1(i)-yd0 
               bspl=h*k1(i)-ydiff 
               cont(j)=y(i) 
               cont(nrd+j)=ydiff 
               cont(2*nrd+j)=bspl 
               cont(3*nrd+j)=-h*k2(i)+ydiff-bspl 
   43       continue 
         end if 
         do  i=1,n 
            k1(i)=k2(i) 
            y(i)=y1(i) 
         end do 
         xold=x 
         x=xph 
         if (iout.ne.0) then 
            hout=h 
            call solout(naccpt+1,xold,x,y,n,cont,icomp,nrd,             &
     &           irtrn)                                       
            if (irtrn.lt.0) goto 79 
         end if 
! ------- normal exit                                                   
         if (last) then 
            h=hnew 
            idid=1 
            return 
         end if 
         if(abs(hnew).gt.hmax)hnew=posneg*hmax 
         if(reject)hnew=posneg*min(abs(hnew),abs(h)) 
         reject=.false. 
      else 
! --- step is rejected                                                  
         hnew=h/min(facc1,fac11/safe) 
         reject=.true. 
         if(naccpt.ge.1)nrejct=nrejct+1 
         last=.false. 
      end if 
      h=hnew 
      goto 1 
! --- fail exit                                                         
   76 continue 
      idid=-4 
      return 
   77 continue 
      if (iprint.gt.0) write(iprint,979)x 
      if (iprint.gt.0) write(iprint,*)' step size t0o small, h=',h 
      idid=-3 
      return 
   78 continue 
      if (iprint.gt.0) write(iprint,979)x 
      if (iprint.gt.0) write(iprint,*)                                  &
     &     ' more than nmax =',nmax,'steps are needed'                  
      idid=-2 
      return 
   79 continue 
      if (iprint.gt.0) write(iprint,979)x 
  979 format(' exit of dopri5 at x=',e18.4) 
      idid=irtrn
      return 
      END                                           
!                                                                       
      function hinit(n,fcn,x,y,xend,posneg,f0,f1,y1,iord,               &
     &     hmax,atol,rtol,itol)                               
! ----------------------------------------------------------            
! ----  computation of an initial step size guess                       
! ----------------------------------------------------------            
      implicit double precision (a-h,o-z) 
      dimension y(n),y1(n),f0(n),f1(n),atol(*),rtol(*) 
! ---- compute a first guess for explicit euler as                      
! ----   h = 0.01 * norm (y0) / norm (f0)                               
! ---- the increment for explicit euler is small                        
! ---- compared to the solution                                         
      dnf=0.0d0 
      dny=0.0d0 
      atoli=atol(1) 
      rtoli=rtol(1) 
      if (itol.eq.0) then 
         do  i=1,n 
            sk=atoli+rtoli*abs(y(i)) 
            dnf=dnf+(f0(i)/sk)**2 
            dny=dny+(y(i)/sk)**2 
         enddo 
      else 
         do i=1,n 
            sk=atol(i)+rtol(i)*abs(y(i)) 
            dnf=dnf+(f0(i)/sk)**2 
            dny=dny+(y(i)/sk)**2 
         end do 
      end if 
      if (dnf.le.1.d-10.or.dny.le.1.d-10) then 
         h=1.0d-6 
      else 
         h=sqrt(dny/dnf)*0.01d0 
      end if 
      h=min(h,hmax) 
      h=sign(h,posneg) 
! ---- perform an explicit euler step                                   
      do  i=1,n 
         y1(i)=y(i)+h*f0(i) 
      enddo 
      call fcn(n,x+h,y1,f1) 
! ---- estimate the second derivative of the solution                   
      der2=0.0d0 
      if (itol.eq.0) then 
         do  i=1,n 
            sk=atoli+rtoli*abs(y(i)) 
            der2=der2+((f1(i)-f0(i))/sk)**2 
            enddo 
         else 
            do  i=1,n 
               sk=atol(i)+rtol(i)*abs(y(i)) 
             der2=der2+((f1(i)-f0(i))/sk)**2 
               end do 
            end if 
            der2=sqrt(der2)/h 
! ---- step size is computed such that                                  
! ----  h**iord * max ( norm (f0), norm (der2)) = 0.01                  
            der12=max(abs(der2),sqrt(dnf)) 
            if (der12.le.1.d-15) then 
               h1=max(1.0d-6,abs(h)*1.0d-3) 
            else 
               h1=(0.01d0/der12)**(1.d0/iord) 
            end if 
            h=min(100*abs(h),h1,hmax) 
            hinit=sign(h,posneg) 
            return 
      END                                           
!                                                                       
      function contd5(ii,x,con,icomp,nd) 
! ----------------------------------------------------------            
!     this function can be used for continuous output in connection     
!     with the output-subroutine for dopri5. it provides an             
!     approximation to the ii-th component of the solution at x.        
! ----------------------------------------------------------            
      implicit double precision (a-h,o-z) 
      dimension con(5*nd),icomp(nd) 
      common /condo5/xold,h 
! ----- compute place of ii-th component                                
      i=0
      do 5 j=1,nd 
         if (icomp(j).eq.ii) i=j 
5        continue
         

      if (i.eq.0) then 
         write (6,*) ' no dense output available for comp.',ii 
         return 
      end if
      theta=(x-xold)/h 
      theta1=1.d0-theta
      contd5=con(i)+theta*(con(nd+i)+theta1*(con(2*nd+i)+theta*         &
           &     (con(3*nd+i)+theta1*con(4*nd+i))))
      return 
      END                                           
!                                                                       
      subroutine cdopri(c2,c3,c4,c5,e1,e3,e4,e5,e6,e7,                  &
     &                    a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,      &
     &                    a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,      &
     &                    d1,d3,d4,d5,d6,d7)                            
! ----------------------------------------------------------            
!     runge-kutta coefficients of dormand and prince (1980)             
! ----------------------------------------------------------            
      implicit double precision (a-h,o-z) 
      c2=0.2d0 
      c3=0.3d0 
      c4=0.8d0 
      c5=8.d0/9.d0 
      a21=0.2d0 
      a31=3.d0/40.d0 
      a32=9.d0/40.d0 
      a41=44.d0/45.d0 
      a42=-56.d0/15.d0 
      a43=32.d0/9.d0 
      a51=19372.d0/6561.d0 
      a52=-25360.d0/2187.d0 
      a53=64448.d0/6561.d0 
      a54=-212.d0/729.d0 
      a61=9017.d0/3168.d0 
      a62=-355.d0/33.d0 
      a63=46732.d0/5247.d0 
      a64=49.d0/176.d0 
      a65=-5103.d0/18656.d0 
      a71=35.d0/384.d0 
      a73=500.d0/1113.d0 
      a74=125.d0/192.d0 
      a75=-2187.d0/6784.d0 
      a76=11.d0/84.d0 
      e1=71.d0/57600.d0 
      e3=-71.d0/16695.d0 
      e4=71.d0/1920.d0 
      e5=-17253.d0/339200.d0 
      e6=22.d0/525.d0 
      e7=-1.d0/40.d0 
! ---- dense output of shampine (1986)                                  
      d1=-12715105075.d0/11282082432.d0 
      d3=87487479700.d0/32700410799.d0 
      d4=-10690763975.d0/1880347072.d0 
      d5=701980252875.d0/199316789632.d0 
      d6=-1453857185.d0/822651844.d0 
      d7=69997945.d0/29380423.d0 
      return 
      END                                           
