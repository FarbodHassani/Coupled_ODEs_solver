implicit none
integer ( kind = 4 ), parameter :: n = 1024*2*2*2*2
  real (kind=8) A
  real (kind=8) B
  real (kind=8) granularity
  integer (kind=4) which
  common/paras/A,B,granularity,which,uuprim,nonlinear
  
  integer*4 shape
  integer*4 extrema
  integer*4 both
  integer*4 derivs
  parameter(shape=0)
  parameter(extrema=1)
  parameter(both=2)
  parameter(derivs=3)
  integer*4 nonlinear
  real (kind=8) , parameter :: pi=4D0*datan(1D0)
  integer*4 uuprim
  
  real (kind=8), parameter ::  maxlimit=1000;
