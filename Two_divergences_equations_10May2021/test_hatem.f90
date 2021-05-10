program test_hatem
  !! gfortran -O6 -ffree-line-length-0 test_hatem.f90 fcn_hatem.f90 dopri5_no_rpar.f90 
  !! ifort -xHost -O3 -w  -ipo test_hatem.f90 fcn_hatem.f90 dopri5_no_rpar.f90 -o hatem
  !! /opt/intel/oneapi/advisor/2021.1.1/bin64/advisor-gui
  !! 
  !!hatem 0.001 0.00003 50 1 2 0 | xmgrace -   

  ! hatem 0.001 0.0001 50 1 2 0 | xmgrace - -world 0 0 4 10

  ! hatem 0.0001 0.001 50 1 2 0 | xmgrace - -world 0 0 4 10  
  !! hatem 0.00001 0.001 50 1 2 0 | xmgrace - -world 0 0 4 10 no divergence doesnt hit bndry
  !! hatem 0.001 0.001 50 1 2 0 0 | xmgrace - -world 0 0 4 10  divergence sideways
  !! hatem 0.001 0.001 50 1 2 0 0 3 | xmgrace - -world 0 0 4 10  central divergence



  !! hatem 0.006 0.0 50 1 2 1 0 3 | xmgrace - -world 0 -10 6 10
  !! divergence changes when b ~0.03
  !! NB this is u'(x,0)\ne 0

  !!  hatem 0.05 0.0 50 0.1 2 0 0 3 | xmgrace - -world 0 -10 6 10
  !! this is u(x,0) ne0

  

!!    y(i+n*beginning)=-exp(-30*(2*pi*(i-n/2)/n)**2)

!!  hatem 0.05 0.15 50 0.2 2 0 0 3 | xmgrace - -world 0 -10 6 10
  !! hatem 0.05 0.2 50 0.2 2 0 0 3 | xmgrace - -world 0 -10 6 10

  !! best has transition between flat and divergence.
!! hatem 0.0008 0.05 50 0.5 2 0 0 3 | xmgrace - -world 0 -10 6 10
!! exp(-30) 16000 pts
  !! transition between
  
!! hatem 0.01 0.05 50 0.5 2 0 0 3 | xmgrace - -world 0 -10 6 10
  !! and
!!   hatem 0.01 0.05 50 0.5 2 0 0 3 | xmgrace - -world 0 -10 6 10
  
  !! n is the number of variables for u and v
  !! nn =2*n is the number of variables in the integration

  include 'fcn.h'
  integer*4 nn
  !! local stuff
  integer*4 i
  integer*4 ndgl,nrdens,lwork,liwork,iwork
  real*8 work,y
  parameter(nn=2*n)
  parameter (ndgl=nn,nrdens=nn)   !! put nrdense=n if iout=2 
  parameter (lwork=8*ndgl+5*nrdens+21,liwork=nrdens+21)

  dimension y(ndgl),work(lwork),iwork(liwork)
  external fcn,solout
  real*8 atol(1),rtol(1),h,x,xend
  integer*4 iout,idid,itol
  character, parameter :: sq = "'"
  real*8 precision
  parameter(precision=1D-11)
  real*8 ux(n),uxx(n)







  !! n is the number of variables for u and v
  !! nn =2*n is the number of variables in the integration





  integer*4 numarg
  real*8 myarg
  external myarg
  integer*4 beginning



  numarg = iargc ( )
  if(numarg.ne.8)then
     print *,"need 8 args: 1:A, 2:B, 3:xend, 4:granularity"
     print *,"5: want to print: u=0,uprime=1,both=2"
     print *,"6: 0 means start with u ne 0, 1 means with uprime ne 0"
     print *,"7: 0 means usual calculation, 1 only nonlinear contribution"
     print *,"8: plot shape=0, extrema=1, u,ux,uxx=2"
     stop
  endif
  A =myarg(1)
  B=myarg(2)
  xend=myarg(3)
  granularity=myarg(4)
  uuprim=myarg(5)
  beginning=myarg(6)  ! 0=u, 1=u'
  nonlinear=myarg(7)
  which=myarg(8)
  if(uuprim.gt.2)then 
     print *,"wrong uuprime parameter, must be 0 for u, 1 for u', 2 for both"
     stop
  endif
  if(beginning.gt.1)then 
     print *,"wrong start parameter, must be 0(for u) or 1 (for uprim"
     stop
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(which==derivs)then
     uuprim=-1
     write(6,'(1h@,"version 50122")')
     write(6,'(1h@,"legend char size 0.460000")')
     write(6,'(1h@,"legend 1.6, 0.9")')
     write(6,'(1h@,"xaxis label ",1h","u      and       ux       and uxx",1h")')
     write(6,'(1h@,"    yaxis  tick major 1")')
     write(6,'(1h@," yaxis  tick minor ticks 0")')
     if(nonlinear==0)then
        write(6,&
             &'(1h@,"title ",1h","A=",E15.6," B=",E15.6," n=",i6," prec=",D8.1,1h")')&
             &A,B,n,precision
     else
        write(6,&
             &'(1h@,"title ",1h","nonlinear part only:  A=",E15.6," B=",E15.6," n=",i6," prec=",D8.1,1h")')&
             &A,B,n,precision
     endif
     write(6,'(1h@,"view 0.100000, 0.150000, 1.75000, 0.850000")')
  endif


  if(uuprim==both)then
     write(6,'(1h@,"version 50122")')
     write(6,'(1h@,"legend char size 0.460000")')
     write(6,'(1h@,"legend 1.3, 0.95")')
     write(6,'(1h@,"xaxis label ",1h","u      and       u_t",1h")')
     which=both
     write(6,'(1h@,"    yaxis  tick major 1")')
     write(6,'(1h@," yaxis  tick minor ticks 0")')
     write(6,&
          &'(1h@,"title ",1h","A=",F10.6," B=",F10.6," n=",i6," prec=",D8.1,1h")')&
          &A,B,n,precision
  endif
  if(uuprim==0)then
     write(6,&
          &'(1h@,"title ",1h","A=",F10.6," B=",F10.6," n=",i6," prec=",D8.1,"u",1h")')&
          &A,B,n,precision
  endif
  if(uuprim==1)then
     write(6,&
          &'(1h@,"title ",1h","A=",F10.6," B=",F10.6," n=",i6," prec=",D8.1," u",a1,1h")')&
          & A,B,n,precision,sq
  endif


  !! take u(x,0)= cos
  !!      u'(x,0)=0
  y=0d00
  do i=1,n
!!    y(i)=1-cos(4*2d0*pi*(i-1)/n)
     !     y(n+n+1-i)=y(i+n)
     y(i+n*beginning)=-exp(-30*(2*pi*(i-n/2)/n)**2)
!     y(i)=0.001*(2*pi*(i-n/2)/n)**2
     y(i+(1-beginning)*n)=0d00

  enddo
  call derivatives(y,ux,uxx)
  do i=1,n
!     y(i+n)=ux(i)

  enddo
 
 
  ! if(which==shape)then
  !    write(6,'(1h@"s0 legend ",1h",F6.3,1h")')0.00
  !    call print_r(y(1+n*uuprim))
  ! endif

  !! prepare dopri5 
  h=0.3d0 
  rtol(1)=precision
  atol(1)=precision
  itol=0    ! rtol and atol are scalars
  iout=2    ! call solout after each step =1, 
  !! iout=2 allows for the interpolation

  ! ----------- default values for dopri5 ------
  do  i=1,20 
     iwork(i)=0 
     work(i)=0.d0
  end do

  !! set start and end
  x=0



  iwork(1)=1000000   !! maximum number of integration steps
  iwork(5)=nn         !! dense output for all variables
  iwork(3)=-1        !! no printing when -1, 0=printing
  !! call integrator
  call dopri5(nn,fcn,x,y,xend,&
       rtol,atol,itol,&
       solout,iout,&
       work,lwork,iwork,liwork,idid)
  if(idid==-3)then
     write(6,'(1h@,"subtitle ",1h","stepsize too small at time t=",f10.6,1h")')x
  endif


end program test_hatem

subroutine solout(nr,xold,x,y,ntotal,con,icomp,nd,irtrn) 
  include'fcn.h'
  integer*4 nr,ntotal,icomp,nd,irtrn
  real*8 xold,x,y,con
  real*8 l2norm
  integer*4 i,it
  integer*4 npts
  real*8 contd5
  real*8 lastdone
  save lastdone
  external contd5
  dimension y(ntotal),con(5*nd),icomp(nd)
  real*8 y1(n),y2(n),current(2*n)
  real*8 timegranularity
  real*8 maxabs
  real*8 ux(n),uxx(n),u00(n)
  integer*4 set
  character(10) leftstr
  save set
  timegranularity=granularity
  irtrn=0
  maxabs=0
     if(nr==1)then
     lastdone=-1
     set=0
 
     return
   endif
  maxabs=maxval(dabs(y))
  if(maxabs>maxlimit)then
     irtrn=-1
    write(6,'(1h@,"subtitle ",1h","diverged at time t=",e15.4,", value=",e15.4,">",E15.4,1h")')x,maxabs,maxlimit
     write(0,*)'diverged at t=',x
  
    return
 endif
 !  write(6,'("#diverged",3e15.7)')x,maxval(dabs(y)),maxabs
  it=x/timegranularity
  !! do we reach the boundary?
  if(abs(y(1))>1.0D-4)then
         irtrn=-1
    write(6,'(1h@,"subtitle ",1h","touched boundary at  t=",e15.4,1h")')x
     write(0,*)'touched boundary at t=',x
 endif 

  
  if(it==lastdone)return
  lastdone=it
  if(which==shape .or. which==both .or.which==derivs)then
     write(leftstr,'(i6)')set
     leftstr=adjustl(leftstr)
     write(6,'(1h@"s",a6," legend ",1h",F6.3,1h")')leftstr,it*timegranularity
     set=set+1
  endif
  do i=1,2*n
     current(i)=contd5(i,it*timegranularity,con,icomp,nd)
  enddo
    call u0(it*timegranularity,u00)
    do i=1,n
       current(i)=current(i)-nonlinear*u00(i)
    enddo

  !! to be done better
  if(which==derivs)then
  call derivatives(current,ux,uxx)
  call print3(current,ux,uxx)
  endif
  if(which==shape)then
     call print_r(current(1+n*uuprim))
  endif
  if(which==both)then
     call print_both(current)
  endif
  if(which==extrema)then
     do i=1,n
        y1(i)=current(i+n*uuprim)
     enddo
     write(6,*)it*timegranularity,minval(y1),maxval(y1)
  endif


end subroutine solout
real*8 function myarg(n)
  implicit none
  integer*4 n
  character*60 arg
  real*8 the_result
  call getarg(n,arg)
  read(arg,*)the_result
  myarg=the_result
end function myarg
subroutine print_both(r)
    include 'fcn.h'
  real*8 r(2*n)
  integer*4 i
  do i =1,2*n
     print *,2D00*i/n,r(i)
  enddo
  print *
end subroutine print_both
subroutine print3(u,ux,uxx)
  include 'fcn.h'
  real*8 u(n),ux(n),uxx(n)
  real*8 all3(3*n)

  integer*4 i
  do i=1,n
     all3(i)=u(i)
     all3(i+n)=ux(i)
     all3(i+2*n)=uxx(i)
  enddo
  do i =1,3*n
     print *,2d0*i/n,all3(i)
  enddo
  print *
end subroutine print3
subroutine u0(t,u)
  include 'fcn.h'
  real*8 t
  real*8 u(n)
  real*8 uini
  external uini
  integer*8 i
  real*8 dx,gamma
  dx=2*pi/n
  gamma=sqrt(B)
  do i=1,n
     u(i)=0.5d00*(uini(dx*i-gamma*t)+uini(dx*i+gamma*t))
  enddo
  end subroutine u0
  real*8 function uini(x)
    include 'fcn.h'
    real*8 x
    uini=-exp(-20*(x-pi)**2)
     return
  end function uini
  
