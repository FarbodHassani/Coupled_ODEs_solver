subroutine fcn(nn,x,y,f)
  ! am not using nn  but a global n
  ! this is needed to compute the tables only once
  !! y and f are in x-space
  include 'fcn.h'
  integer*4 nn
  real*8 x,y(n+n),f(n+n)
  logical,save :: firstcall=.true.

  !! local variables
  integer*4 i
  real*8 ux(n),uxx(n)
  real*8 u(n),v(n)

  !! make u and v=u'
  do i=1,n
     u(i)=y(i)
  enddo
  !! dot u=v
  do  i=1,n
     f(i)=y(i+n)
  enddo
  !! now the rhs

  call derivatives(y,uxx,ux)
  do i=1,n
     f(n+i)=A*ux(i)**2+B*uxx(i)
  enddo
end subroutine fcn

subroutine print_r(r)
  include 'fcn.h'
  real*8 r(n)
  integer*4 i
  do i =1,n
     print *,2*pi*i/n,r(i)
  enddo
  print *
end subroutine print_r

subroutine derivatives(rin0,rout,rout2)
  include 'fcn.h'
  real*8  rin(-1:n+2),rout(n),rout2(n),rin0(n)
  real*8 am2,am1,a0,ap1,ap2
  integer*4 i
  real*8 dx
  dx=2d0*pi/n
  do i=1,n
     rin(i)=rin0(i)
  enddo
  !extend periodically
  rin(0)=rin0(n)
  rin(-1)=rin0(n-1)
  rin(n+1)=rin0(1)
  rin(n+2)=rin0(2)

  do i=1,n
     am2=rin(i-2)
     am1=rin(i-1)
     a0=rin(i)
     ap1=rin(i+1)
     ap2=rin(i+2)
     rout2(i)=(1d0*am2-8d0*am1+8d0*ap1-1d0*ap2)/(12d00*dx)
     rout(i)=(-1D0*am2+16d0*am1-30d0*a0+16d0*ap1-1d0*ap2)/(12d00*dx**2)
  enddo
end subroutine derivatives
