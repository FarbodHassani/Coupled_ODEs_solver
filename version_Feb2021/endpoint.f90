program  cosmology
  !!gfortran -O6 -ffree-line-length-0 endpoint dopri5_no_rpar.f90 dopri5copy_no_rpar.f90 -o endpoint
  !! ifort -xHost -O3 -w -qopt-report -prof-use -ipo endpoint.f90 dopri5_no_rpar.f90 dopri5copy_no_rpar.f90  -o endpoint

  !! usage: plotendpoints.pl >! file.txt;
  !! full_prepare_for_matlab0203ampl1minuscs2.pl
  !! or: full_prepare_for_matlab0203ampl.pl*
  !!depending whether you are interested in 1-cs^2 or cs^2

  
  !! usage 
  !! we need 2 distinct copies of dopri5 because we use dopri5copy within dopri5
  !! (Thanks to Ernst Hairer for helping the debug)  
  !! ifort runs under bash
  !! seems about 25% faster

  !! we are using here a modified version of dopri5
  !! in which rpar and ipar are omitted, since we do not use them
  !! the gain in time is practically negligible
  !! but the program becomes a little less unreadable
  include 'physics2.h'
  !! variables for dopri5
  integer*4 ndgl,nrdens,n,lwork,liwork,iwork
  real*8 y,work,cosmofcn
  parameter(n=1+4*(degree+1))
  parameter (ndgl=n,nrdens=n)   !! put nrdense=n if iout=2 
  parameter (lwork=8*ndgl+5*nrdens+21,liwork=nrdens+21) 
  dimension y(ndgl),work(lwork),iwork(liwork)
  external cosmofcn,cosmosolout
  real*8 atol(1),rtol(1),h,x,xend
  integer*4 i,iout,idid,itol

  !! variables for the calculation
  integer*4 paramw,paramc,nw,nc
  real*8 cmin,cmax,wmin,wmax
  real*8 z2tau
  external z2tau
  real*8 tauini,tauend

  real*8 exp
  real*8 psi2_i,cs2_i,w_i  ! the initial values
  character*60 arg1,arg2,arg3,arg4,arg
  integer*4 numarg
  nw=0
  nc=0

  wmin=-1.5d00
  wmax=0.5d00

  cmin=0d00
  cmax=1d00


  !  write(0,*)'psi=',psi+1,'psip=',psip+1,'pi=',pi+1,'pip=',pip+1

  numarg = iargc ( )
  if(numarg.ne.3)then
     print *,"need 3 args: psi2, cs2,w"
     stop
  endif
  call getarg ( 1, arg1 )
  read(arg1,*)psi2_i
  call getarg ( 2, arg2 )
  read(arg2,*)cs2_i
  call getarg ( 3, arg3 )
  read(arg3,*)w_i

!  write(6,*)'psi2=',psi2_i
!  write(6,*)'cs2=',cs2_i
!  write(6,*)'w=',w_i
!    write(6,*)'1-cs2=',1D0-cs2_i,cs2

  kvalue=[0d0,0d0,psi2_i,0d0,0d0]



  do paramw=0,nw   ! loop over w
     !     w=-1.5d00+1.5d00*paramw/nw ! to be changed (should depend on paramw)
     w=w_i

     tauini=z2tau(zini)  ! depends on w
     tauend=z2tau(0d00)  !depends on w 
!          write(0,*)tauini,tauend
     do paramc=0,nc   !loop over c_s^2
        !        cs2=(cmax-cmin)*paramc/nc+cmin
        cs2=cs2_i
        !! the initial conditions        
        do i=1,n
           y(i)=0
        enddo
        ! we set psi(tau_i)=kvalue(i)
        !(everything else is 0)
        !! NB, the initial value of a at tauini is 1/(zini+1)
        y(a)=1/(zini+1)
        ! initializing the psi
        do i=0,degree
           y(psi+i)=kvalue(i)
        enddo
        !        y(pip+2)=0.0d00
        !! prepare dopri5 
        h=0.3d0 
        rtol(1)=precision
        atol(1)=precision
        itol=0    ! rtol and atol are scalars
        iout=2    ! call cosmosolout after each step =1, 
        !! iout=2 allows for the interpolation

        ! ----------- default values for dopri5 ------
        do  i=1,20 
           iwork(i)=0 
           work(i)=0.d0
        end do

        !! set start and end
        x=tauini
        xend=tauend
        

        iwork(1)=1000000   !! maximum number of integration steps
        iwork(5)=n         !! dense output for all variables
        iwork(3)=-1        !! no printing when -1, 0=printing
        !! call integrator
        call dopri5(n,cosmofcn,x,y,xend,&
             rtol,atol,itol,&
             cosmosolout,iout,&
             work,lwork,iwork,liwork,idid)

        write(6,'(3e25.17,i4)')cs2,psi2_i,x,idid
        goto 12345
        if(idid==-15)then
           write(0,*)'singularity detected at taustar',x
        endif
        if(idid==-3)then
           write(0,*)'stepsize too small',x
        endif
        if(idid==-17)then
           write(0,*)'l2norm > 10^{20}',x
        endif
12345 continue

     enddo  !cs2

  enddo !w
end program cosmology


subroutine cosmosolout(nr,xold,x,y,n,con,icomp,nd,irtrn) 
  include'physics2.h'
  integer*4 nr,n,icomp,nd,irtrn
  real*8 xold,x,y,con
  real*8 l2norm
  integer*4 i
  real*8 contd5
  external contd5
  dimension y(n),con(5*nd),icomp(nd)
  real*8 z2tau,tau
  external z2tau
  real*8 interpolate(21)
  real*8 nextz ! the z start at high and descend to 0
  real*8 taustar,taustar2, ampl
  real*8 x1,x2,x3,f1,f2,f3,newtaustar
  external newtaustar
  save nextz
  if(nr==1)then
     nextz=zini
  endif
  irtrn=1
  l2norm=0
  do i=0,4
     l2norm=l2norm+y(pi+i)**2+y(pip+i)**2
  enddo
  if(l2norm>10D40**2)then
     irtrn=-17
     return
  endif
  !! for intermediate values use contd5(i,time,con,icomp,nd)
  !! equispaced at integer values of z (not tau)
!!  do while(1d00/y(a)-1d00 < nextz)
     tau=z2tau(nextz)
     do i=1,nd
        interpolate(i)=contd5(i,tau,con,icomp,nd)
     enddo
!!     write(6,'(22e15.7)')x,y
!!     nextz=nextz-1
!!  enddo
  if(nr>1)then
!     call   ampltaustar(x,xold,y(pi+2),contd5(pi+2,xold,con,icomp,nd),taustar,ampl)
!     call   ampltaustar2(x,xold,y(pi+2),contd5(pi+2,xold,con,icomp,nd),taustar2,ampl)
     x1=x
     x2=x-0.05d00
     x3=x+0.05d00
     f1=contd5(pip+2,x1,con,icomp,nd)
     f2=contd5(pip+2,x2,con,icomp,nd)
     f3=contd5(pip+2,x3,con,icomp,nd)
     
     taustar=newtaustar(x1,x2,x3,f1,f2,f3)
     !!    if(taustar<10*x.and. taustar>0)  print *,x,abs(taustar+0.00001),ampl,y(pi+2)
!!     write(6,'(25e15.7)')x,y,taustar-x,taustar2,l2norm
     if(abs(taustar-x)<0.1D00 .and. abs(y(pip+2))>5)then
        irtrn=-15
     endif
  endif
  return 
end subroutine cosmosolout

subroutine cosmofcn(n,x,y,f)
  include 'physics2.h'
  integer*4 ::  n
  real*8 :: y,x
  real*8  :: f
  dimension y(n),f(n)
  !! the coefficents (some do not depend on x, nor y
  !! optimizer will move them to the relevant place in program
  !! the labels are pi for \pi, psi for  \psi
  !! and attached d means derivative with respect to (conformal) time
  !! p derivative w.r.t. x
  real*8 ell_pi
  real*8 ell_pid
  real*8 ell_pipp
  real*8 ell_psi
  real*8 ell_psid
  real*8 nu_pid_pipp
  real*8 nu_pi_pipp
  real*8 nu_pip_pip
  real*8 nu_pip_pipd
  real*8 nu_pip_pip_pipp
  real*8 nu_psi_pipp
  real*8 nu_psip_pip
  real*8 hx,hprimex
  real*8 hprime,psirhs
  external hprime,psirhs
  
  integer*4 i
!!! these definitions should make the formulas for f(...) more readable.
  !! perhaps they give a factor 3 of speedup
  
  real*8 pi0,pi1,pi2,pi3,pi4
  real*8 pid0,pid1,pid2,pid3,pid4
  real*8 psi0,psi1,psi2,psi3,psi4
  real*8 psid0,psid1,psid2,psid3,psid4

  !! note that the suffix p here really means derivative wrt tau, not x
  !! the x derivatives are captured via the digits 0-4 (which select the
  !! coefficients of the powers of x
  
  pi0=y(pi+0)
  pi1=y(pi+1)
  pi2=y(pi+2)
  pi3=y(pi+3)
  pi4=y(pi+4)
  pid0=y(pip+0)
  pid1=y(pip+1)
  pid2=y(pip+2)
  pid3=y(pip+3)
  pid4=y(pip+4)

  

  psi0=y(psi+0)
  psi1=y(psi+1)
  psi2=y(psi+2)
  psi3=y(psi+3)
  psi4=y(psi+4)
  psid0=y(psip+0)
  psid1=y(psip+1)
  psid2=y(psip+2)
  psid3=y(psip+3)
  psid4=y(psip+4)



  !! the equation for a
  !!  this must be computed first, since we will need f(a) below

  ! a'=h0*sqrt(Odarke*a**(1-3*w)+Omatter*a+Oradiation)
  f(a)=h0*sqrt(Odarke*y(a)**(1-3*w)+Omatter*y(a)+Oradiation)

  !! now that we have a and a' we can compute H and H' (at x)
  !! we write hx (because h is the step size in dopri5) 

  !! Note that H(x) is simply a'(x)/a(x) i.e., f(a)/y(a)
  !! while       H'(x) = Q(a(x))
  !! Q is computed in function hprime
  !! so we dont need an additional integration for H or H'

  hx=f(a)/y(a)
  hprimex=hprime(y(a))
 ! hx=2/x
 ! hprimex=-2/x**2
  
  !! the coefficients which appear in the u_i equations
  !! since x is fixed in cosmofcn, 
  !! but of course they do through the current values of hx,hprimex

  !!the ell
  ell_pid  = (1-3*w)*hx                         ! depends on x
  ell_pi   =(1-3*cs2)*hprimex+3*(cs2-w)*hx**2   ! depends on x
  ell_psi  = 3*(w-cs2)*hx                       ! depends on x
  ell_psid =-(1+3*cs2)
  ell_pipp = -cs2


  !! same for the nu
  nu_pip_pip       = -hx*(5*cs2+3*w-2)/2        ! depends on x
  nu_pip_pipd     = 2*(1-cs2)
  nu_pi_pipp     = -hx*(cs2-1+3*cs2*(1+w))      ! depends on x
  nu_pid_pipp   = -(cs2-1)
  nu_psi_pipp    = (cs2-1)
  nu_psip_pip    = (cs2-1)
  nu_pip_pip_pipp= 3*(cs2-1)/2


  !the equation for psi_0 to psi_4
  do i=0,degree
     f(psi+i)=y(psip+i)                  ! psi_i
  enddo
  do i=0,degree
     f(psip+i)=psirhs(n,y,i,hx,hprimex)   ! psi_i'
  enddo
 

  

  !! now we have all the psi,psi' which appear in the equs for the u_i 
  !! the equations for the u_i follow

  !! the second derivatives u_i'=v_i , i=0,...,4
  do i=0,degree
     f(pi+i)=y(pip+i)
  enddo

  !! these equations are generated from the mathematica output
  !! the first derivative v_i'=F(u,v,...), i=0,...4
  f(pip+0)=-(ell_psi*psi0+ell_psid*psid0+ell_pid*pid0&
       +2*pi2*ell_pipp+ell_pi*pi0-2*pi2*nu_pid_pipp*pid0&
       +pi1**2*(-nu_pip_pip)-2*pi0*pi2*nu_pi_pipp&
       -pi1*nu_pip_pipd*pid1-2*pi2*pi1**2*nu_pip_pip_pipp&
       -2*psi0*pi2*nu_psi_pipp-psi1*pi1*nu_psip_pip)

  f(pip+1)=-(ell_psi*psi1+ell_psid*psid1+ell_pid*pid1+6*pi3*ell_pipp&
       +ell_pi*pi1-6*pi3*nu_pid_pipp*pid0-2*pi2*nu_pid_pipp*pid1&
       -4*pi2*pi1*nu_pip_pip-2*pi2*pi1*nu_pi_pipp&
       -6*pi0*pi3*nu_pi_pipp-2*pi1*nu_pip_pipd*pid2&
       -2*pi2*nu_pip_pipd*pid1-6*pi3*pi1**2*nu_pip_pip_pipp&
       -8*pi2**2*pi1*nu_pip_pip_pipp&
       -6*psi0*pi3*nu_psi_pipp-2*psi1*pi2*nu_psip_pip&
       -2*psi1*pi2*nu_psi_pipp-2*psi2*pi1*nu_psip_pip)

  f(pip+2)=-(ell_psi*psi2+ell_psid*psid2+ell_pid*pid2+12*pi4*ell_pipp&
       +ell_pi*pi2-2*pi2*nu_pid_pipp*pid2&
       -12*pi4*nu_pid_pipp*pid0-6*pi3*nu_pid_pipp*pid1&
       -4*pi2**2*nu_pip_pip-6*pi1*pi3*nu_pip_pip&
       -2*pi2**2*nu_pi_pipp-6*pi1*pi3*nu_pi_pipp&
       -12*pi0*pi4*nu_pi_pipp-4*pi2*nu_pip_pipd*pid2&
       -3*pi3*nu_pip_pipd*pid1-3*pi1*nu_pip_pipd*pid3&
       -8*pi2**3*nu_pip_pip_pipp-36*pi1*pi3*pi2*nu_pip_pip_pipp&
       -12*pi1**2*pi4*nu_pip_pip_pipp-12*psi0*pi4*nu_psi_pipp&
       -3*psi1*pi3*nu_psip_pip-6*psi1*pi3*nu_psi_pipp&
       -4*psi2*pi2*nu_psip_pip-2*psi2*pi2*nu_psi_pipp&
       -3*psi3*pi1*nu_psip_pip)

  f(pip+3)=-(ell_psi*psi3+ell_psid*psid3+ell_pid*pid3+ell_pi*pi3&
       -2*pi2*nu_pid_pipp*pid3-12*pi4*nu_pid_pipp*pid1&
       -6*pi3*nu_pid_pipp*pid2-12*pi3*pi2*nu_pip_pip&
       -8*pi1*pi4*nu_pip_pip-8*pi3*pi2*nu_pi_pipp&
       -12*pi1*pi4*nu_pi_pipp-6*pi2*nu_pip_pipd*pid3&
       -4*pi4*nu_pip_pipd*pid1-6*pi3*nu_pip_pipd*pid2&
       -4*pi1*nu_pip_pipd*pid4-48*pi3*pi2**2*nu_pip_pip_pipp&
       -64*pi1*pi4*pi2*nu_pip_pip_pipp&
       -36*pi1*pi3**2*nu_pip_pip_pipp&
       -4*psi1*pi4*nu_psip_pip-12*psi1*pi4*nu_psi_pipp&
       -6*psi2*pi3*nu_psip_pip-6*psi2*pi3*nu_psi_pipp&
       -6*psi3*pi2*nu_psip_pip-2*psi3*pi2*nu_psi_pipp&
       -4*psi4*pi1*nu_psip_pip)

  f(pip+4)=-(ell_psi*psi4+ell_psid*psid4+ell_pid*pid4+ell_pi*pi4&
       -2*pi2*nu_pid_pipp*pid4-12*pi4*nu_pid_pipp*pid2&
       -6*pi3*nu_pid_pipp*pid3-16*pi4*pi2*nu_pip_pip&
       -9*pi3**2*nu_pip_pip-14*pi4*pi2*nu_pi_pipp&
       -6*pi3**2*nu_pi_pipp-8*pi2*nu_pip_pipd*pid4&
       -8*pi4*nu_pip_pipd*pid2-9*pi3*nu_pip_pipd*pid3&
       -80*pi4*pi2**2*nu_pip_pip_pipp&
       -90*pi3**2*pi2*nu_pip_pip_pipp&
       -120*pi1*pi3*pi4*nu_pip_pip_pipp&
       -8*psi2*pi4*nu_psip_pip-12*psi2*pi4*nu_psi_pipp&
       -9*psi3*pi3*nu_psip_pip-6*psi3*pi3*nu_psi_pipp&
       -8*psi4*pi2*nu_psip_pip-2*psi4*pi2*nu_psi_pipp)

end subroutine cosmofcn


real*8 function z2tau(z)  !! computes tau as function of z
  !! NB. This needs an integration 
  include 'physics2.h'

  ! note that tau(z_*=infty,w) =0
  ! therefore we can write the diffential equation
  ! tau'(z) =-H_0^{-1} 1/sqrt(...) *

  ! numerically, wew should not integrate to infinity, but 
  ! we can bound the contribution int_Z^\infty by
  ! using only the term 1/\sqrt{\Omega_m (z)^3}
  !  with Omega_m=0.31 and H0^{-1}=15
  ! < 51/sqrt(Z). so, if we want to know tau (which is at least 0.01 or so)
  ! we take 51/sqrt(Z)=0.00001
  ! which is about Z= 10^{16}
  ! 


  real*8 z

  !! dopri5 stuff

  integer*4 ndgl,nrdens,n,lwork,liwork,iwork
  real*8 y,work
  parameter(n=1)
  parameter (ndgl=n,nrdens=n)   !! put nrdense=n if iout=2 
  parameter (lwork=8*ndgl+5*nrdens+21,liwork=nrdens+21) 
  dimension y(ndgl),work(lwork),iwork(liwork)
  real*8 h,atol(1),rtol(1),x,xend
  integer*4 itol,iout,idid
  external taufcn,tausolout

  !! other variables
  integer*4 i




  y=0D00
  x=z ! this is the initial point of the integration 
  xend=1.0D16 ! this is a sufficiently large upper bound to guarantee
  ! tau to error 0.000001

!!!prepare dopri5
  h=0.3d0 
  rtol(1)=precision
  atol(1)=precision
  itol=0    ! rtol and atol are scalars
  iout=0    ! call tausolout after each step
  !! iout=2 allows for the interpolation
  !! we dont need this so far

  ! ----------- default values ------
  do  i=1,20 
     iwork(i)=0 
     work(i)=0.d0
  end do
  iwork(1)=1000000
  iwork(5)=n     ! in case we need dense output
  iwork(3)=-1   !no printing
  call dopri5copy(n,taufcn,x,y,xend,&
       rtol,atol,itol,&
       tausolout,iout,&
       work,lwork,iwork,liwork,idid)
   if(idid.ne.1)write(0, *)'# failure in z2tau for w=',w,' idid=',idid
  z2tau=y(1)
  return
end function z2tau

subroutine tausolout (nr,xold,x,y,n,con,icomp,nd,irtrn) 
  !! currently not used
  implicit none
  integer*4 nr,n,icomp,nd,irtrn
  real*8 xold,x,y,con
  real*8 contd5copy
  external contd5copy
  dimension y(n),con(5*nd),icomp(nd)
  real*8 f(n)
  irtrn=1
  if(nr.ge.1)then
     call taufcn(n,x,y,f)
     print *,x,y(1),nr
  end if
  !! for intermediate values use contd5copy(i,time,con,icomp,nd)
  return 
end subroutine tausolout


subroutine taufcn(n,x,y,f)
  include 'physics2.h'
  integer*4 n
  real*8 y,f,x
  dimension y(n),f(n)
  f(1)=h0inv/sqrt(Odarke*(1+x)**(3*(1+w))+Omatter*(1+x)**3+Oradiation*(1+x)**4)
end subroutine taufcn
  
real*8 function hprime(aa)
  include 'physics2.h'
  real*8 aa ! this is a(tau)
  real*8 t1,t2,t3
  t1=Odarke*(3*w+1)*aa**(-1-3*w)/2d00
  t2=Oradiation/aa**2
  t3=Omatter/(2*aa)
  hprime= -h0**2*(t1+t2+t3)
end function hprime

real*8 function psirhs(n,y,i,hx,hprimex)
  include 'physics2.h'
  integer*4 i,n
  real*8 y(n),hx,hprimex
  psirhs=-3*hx*y(psip+i)&
      -((2-3*Omatter/(2*Omatter+2*Odarke/y(a)**(3*w)+2*Oradiation/y(a)))*hx**2&
      +hprimex)*y(psi+i)
end function psirhs
 
subroutine ampltaustar(x,xold,u,uold,ampl,taustar)
  include 'physics2.h'
  real*8 x,xold,u,uold,ampl,taustar
  real*8 s
  s=(x-xold)*sqrt(abs(u*uold))
  taustar=(u*x-uold*xold+s)/(u-uold)
  ampl=(u*uold*(x-xold)*(2*s+(x-xold)*(u+uold)))/(u-uold)**2
  
end subroutine ampltaustar
subroutine ampltaustar2(x,xold,u,uold,ampl,taustar)
  include 'physics2.h'
  real*8 x,xold,u,uold,ampl,taustar
  real*8 s
  s=(x-xold)*sqrt(abs(u*uold))
  taustar=(u*x-uold*xold)/(u-uold)
  ampl=(u*uold*(x-xold))/(u-uold)
  
end subroutine ampltaustar2
real*8  function newtaustar(x1,x2,x3,f1,f2,f3)
  real*8 x1,x2,x3,f1,f2,f3
   newtaustar = (((x1 ** 2 - x2 ** 2) * f2 ** 2 + (-x1 ** 2 + x3 ** 2) * f3 &
     ** 2) * f1 ** 2 + f3 ** 2 * f2 ** 2 * (x2 - x3) * (x2 + x3)) / (((&
     2 * x1 - 2 * x2) * f2 ** 2 - 2 * f3 ** 2 * (x1 - x3)) * f1 ** 2 + &
     2 * (x2 - x3) * f3 ** 2 * f2 ** 2)
 end function newtaustar
 
