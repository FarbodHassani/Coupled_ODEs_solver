program test_fcn

  !! /Users/farbod/anaconda2/bin/gfortran -c fftpack5.f90
  !! gfortran -O6 -ffree-line-length-0 test_extrema.f90 fcn_extrema.f90 dopri5_no_rpar.f90 fftpack5.o
  !! ifort -xHost -O3 -w  -ipo test_extrema.f90 fcn_extrema.f90 dopri5_no_rpar.f90 fftpack5.f90
  !! /opt/intel/oneapi/advisor/2021.1.1/bin64/advisor-gui
  !!
  !! a.out 0.6 1 | grep # | sed -e's/#//' | xmgrace -




  !! n is the number of variables for u and v
  !! nn =2*n is the number of variables in the integration

  include 'fcn.h'
  integer*4 nn
  !! local stuff
  real*8 pi
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
  real*8 precision
  parameter(precision=1D-8)




  integer*4 numarg
  real*8 myarg
  external myarg


  numarg = iargc ( )
  if(numarg.ne.3)then
     print *,"need 3 args: A, B, xend"
     stop
  endif
  A =myarg(1)
  B=myarg(2)
  xend=myarg(3)

  write(6,'(1h@,"title ",1h","A=",F5.3," B=",F5.3," n=",i6," prec=",D8.1,1h")')A,B,n,precision

!  write(6,'('@title "A=','xend=',xend,'n=',n,"prec=",precision," "
   print *, 0,0
  print *, 2*3.1415926,0
  print *




  pi=4d00*atan(1d00)
  !! take u(x,0)= cos
  !!      u'(x,0)=0
  do i=1,n
     y(i+n)=-cos(2*pi*(i-1)/n)
!    y(i)=1-exp(-10*(2*pi*(i-n/2)/n)**2)
     y(i)=0d00
  enddo
!  call print_r(y,n)


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



end program test_fcn

subroutine solout(nr,xold,x,y,ntotal,con,icomp,nd,irtrn)
  include'fcn.h'
  integer*4 nr,ntotal,icomp,nd,irtrn
  real*8 xold,x,y,con
  real*8 l2norm
  integer*4 i
  integer*4 npts
  real*8 contd5
  external contd5
  dimension y(ntotal),con(5*nd),icomp(nd)
  real*8 y1(n),y2(n)
  irtrn=0

 ! if(int(nr/10000)*10000.ne.nr)return
  if(nr==1)then
  endif
  do i=1,n
     y1(i)=y(i)
     y2(i)=y(i+n)
  enddo
  write(6,*)x,minval(y1),maxval(y1)
!  if(minval(y1)<-0.0001d00)irtrn=-1
  npts=ntotal/2
!  call print_r(y,npts)
  if(maxval(abs(y))>1D10)then
     write(0,'("diverged",e15.7)')maxval(abs(y))
     irtrn=-1
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
