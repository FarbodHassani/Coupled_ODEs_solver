subroutine fcn(nn,x,y,f)
  ! am not using nn  but a global n
  ! this is needed to compute the tables only once
  !! y and f are in x-space
  include 'fcn.h'
  integer*4 nn
  real*8 x,y(n+n),f(n+n)
  logical,save :: firstcall=.true.

  complex ( kind = 8 ) c(n),cprime2(n),cprime(n)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ),parameter :: inc=1
  integer ( kind = 4 ), parameter :: lenc=n
  integer ( kind = 4 ), parameter :: lensav= 2 * n + int ( log ( real ( n, kind = 8 ) ) ) + 4 
  integer ( kind = 4 ), parameter::  lenwrk=2*n
  real ( kind = 8 )  work(1:lenwrk)
  real ( kind = 8 ) wsave(1:lensav)
  save wsave
  !! local variables
  integer*4 i
  real*8 pi
  real*8 der(n)
  complex ( kind = 8 )u(n),uk(n),ukprime(n),ukprimeprime(n)
  real (kind = 8)uxsquare(n),uxx(n)

  !! we only initialize tables once
  if(firstcall)then
     firstcall=.false.
     call zfft1i ( n, wsave, lensav, ier )
     endif
  !! make u and v=u'
  do i=1,n
     u(i)=y(i)
  enddo
  !! dot u=v
  do  i=1,n
     f(i)=y(i+n)
  enddo
  !! now the rhs
  !! convert u to uk (the fourier transform)
  do i=1,n
     uk(i)=y(i)
  enddo
  call zfft1f ( n, inc, uk, lenc, wsave, lensav, work, lenwrk, ier )
  call test(ier)
  
  !now multiply by I k
  call derivative(uk,ukprime)
!  call print_c(ukprime,n)
!  stop
  !! and transform back
  call zfft1b ( n, inc, ukprime, lenc, wsave, lensav, work, lenwrk, ier )
  call test(ier)
  !! now ukprime is actually in x-space
  !!  square it
  !at this point, ukprime**2 should be real

  do i=1,n
     uxsquare(i)=ukprime(i)*ukprime(i)
  enddo

 
  !! u''  (here we could gain a little bit by using that uk is actually real)
  call derivative2(uk,ukprimeprime)
    !! and transform back
  call zfft1b ( n, inc, ukprimeprime, lenc, wsave, lensav, work, lenwrk, ier )
  !! this should be real
  do i=1,n
     uxx(i)=real(ukprimeprime(i),kind=8)*(1d00+2d0/n)
  enddo
!  call derivative2a(y,der)
  do i=1,n
!     print *,'###',i,uxx(i)/der(i)
  enddo
  
  !! finally, we fill f

  do i=1,n
     f(n+i)=A*uxsquare(i)+B*uxx(i)
  enddo

end subroutine fcn
subroutine derivative(cin,cout)
include 'fcn.h'
  complex ( kind = 8 ) cin(n),cout(n)
  real*8 pi
  integer*4 i
  real*8 q
  q=2d00/n

  
 do i=1,n
   cout(i)=(q*(i-1))*cmplx(-dimag(cin(i)),real(cin(i),kind=8))
 enddo
end subroutine derivative
subroutine derivative2(cin,cout)
  include 'fcn.h'
  complex ( kind = 8 ) cin(n),cout(n)
  real*8 pi
  integer*4 i
  real*8 q
  
  q=2d00/n
 do i=1,n
    !    cout(i)=-((2d0*(i-n/2d0)/n)**2)*cin(i)
        cout(i)=-(((i-1)*q)**2)*cin(i)
 enddo
end subroutine derivative2



subroutine test(ier)
  implicit none
  integer*4 ier
  if(ier.ne.0)then
     print *,ier
     stop
  endif
end subroutine  test
subroutine print_r(r,n)
  implicit none
  integer*4 n
  real*8 r(n)
  integer*4 i
  real*8 pi
  return
  pi=4d00*atan(1D00)
  do i =1,n
     print *,2*pi*i/n,r(i)
  enddo
  print *
end subroutine print_r
subroutine print_c(c,n)
  implicit none
  integer*4 n
  complex ( kind = 8 ) c(n)
  integer*4 i
  real*8 pi

  pi=4d00*atan(1D00)
  do i =1,n
     print *,2*pi*i/n,real(c(i),kind=8)
  enddo
  print *
    do i =1,n
     print *,2*pi*i/n,dimag(c(i))
  enddo
  print *
end subroutine print_c
subroutine derivative2a(rin,rout)
  include 'fcn.h'
  real*8  rin(n),rout(n)
  real*8 aa,bb,c,dx
  integer*4 i
  real*8 pi
  pi=4d0*atan(1d0)
  dx=2*pi/n
  aa=rin(n)
  bb=rin(1)
  c=rin(2)
  rout(1)=4D00 *(aa+c-2*bb)/(2*dx**2)
  aa=rin(n-1)
  bb=rin(n)
  c=rin(1)
  rout(n)=4D00 *(aa+c-2*bb)/(2*dx**2)
!  print *,rout(1),rout(n)
  do i=2,n-1
     aa=rin(i-1)
     bb=rin(i)
     c=rin(i+1)
     rout(i)=4D00 *(aa+c-2*bb)/(2*dx**2)
  enddo
end subroutine derivative2a
