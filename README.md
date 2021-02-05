# Fortran_Coupled_ODEs_solver
Fortran code and perl to solve coupled k-essence ODEs written by Jean-Pierre Eckmann using Ernst Hairer (http://www.unige.ch/~hairer/) solver, this code validated by Farbod Hassani's mathematica notebook. 

To run:
gfortran -O6 -ffree-line-length-0 one_orbit.f90 dopri5_no_rpar.f90 dopri5copy_no_rpar.f90 -o one_orbit 
or
ifort -xHost -O3 -w -qopt-report -prof-use -ipo one_orbit.f90 dopri5_no_rpar.f90 dopri5copy_no_rpar.f90  -o one_orbit
usage : one_orbit psi2 cs2 w | plot_one.pl | xmgrace -

The cosmology is set in physics2.h. The initial and final time is set in the one_orbit.f90
The output is a text file with many columns corresponsing to the pi_0, pi_1, pi_2 ..., psi_0, psi_1 ...pi'_0 ....
To compile we can write:
./one_orbit 1 1 1

To plot we can use the perl file:
Dont forget to make the perls excecutable. 


