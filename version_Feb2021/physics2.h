implicit none
real*8 Oradiation,Omatter,Odarke,h0,h0inv,w,cs2,zini,zfinal,precision
  parameter(Oradiation=9D-5)
  parameter(Omatter=0.31D00)
 !! parameter(Omatter=1D0-Oradiation)     !no dark energy
  parameter(Odarke=1d00-Oradiation-Omatter) ! the sum should be exactly 1
  parameter(h0=0.0691023D00)
  parameter(h0inv=1d00/h0)
  parameter(zini=1000D00)
!! parameter(zfinal=-0.9999D00)
   parameter(zfinal=0D00)
  parameter(precision=0.3D-12)
  
  !! these are the initial k values for the
  !! coefficients of of x^j of psi 
  real*8 kvalue(0:4)


  
  common /parameters/w,cs2,kvalue

  !! The vector y has 1+4*(degree+1) components
  !! We number them symbolically to make the program less unreadable
  !! the indices into y
  integer*4 a,pi,pip,psi,psip
  integer*4 degree
  parameter(degree=4)  !!
  !! please note that the program cannot be used for any degree othe that 4
  !! without substantial changes in f(u+xx)=... (in routine cosmofcn)
  parameter(a=1)
  parameter(psi=2) !there are psi_0,...,psi_4
  parameter(psip=2+degree+1)
  parameter(pi=2+2*(degree+1))  !there are u_0,...,u_4
  parameter(pip=2+3*(degree+1))
  
  

!! saving data for output
  real*8 lastx,lasttaustar,lasty,lastampl
  common /output/ lastx,lasttaustar,lasty,lastampl
 
