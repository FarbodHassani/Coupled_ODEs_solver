# Coupled ODEs solver

Fortran code and perl to solve coupled k-essence ODEs written by Jean-Pierre Eckmann using Ernst Hairer (http://www.unige.ch/~hairer/) solver, this code validated by Farbod Hassani's mathematica notebook. The code includes many commments thakns to Jean-Pierre Eckamnn.

## Single orbit
```bash
gfortran -O6 -ffree-line-length-0 one_orbit.f90 dopri5_no_rpar.f90 dopri5copy_no_rpar.f90 -o one_orbit 
```
or
```bash
ifort -xHost -O3 -w -qopt-report -prof-use -ipo one_orbit.f90 dopri5_no_rpar.f90 dopri5copy_no_rpar.f90  -o one_orbit
```
usage :
```bash
one_orbit psi2 cs2 w | plot_one.pl | xmgrace -
```

The cosmology is set in physics2.h. The initial and final time is set in the one_orbit.f90 file.
The output is a text file with many columns corresponsing to the pi_0, pi_1, pi_2 ..., psi_0, psi_1 ...pi'_0 ....
To compile:

```bash
./one_orbit 1 1 1
```

To plot we can use the perl file:
Dont forget to make the perls excecutable. 

```bash
./one_orbit 1 1 -0.7 | ./plot_one.pl | xmgrace -
```

Also to save data:

```bash
./one_orbit 1 1 -0.7 | ./plot_one.pl > data.txt
```
or  with ">!" to rewrite the availbe data file.

```bash
./one_orbit 1 1 -0.7 | ./plot_one.pl >! data.txt
```

In xmgrace -> Plot --> set appearence --> choose data set and in type to see how to represent the data and also check the solout precision/time steps.



## Study for varied quantities

For the grid runs and going over different values of c_s^2 and w:
It's better to use tcsh shell: /bin/tcsh to have access to more commands

First compile the fortran with either of the following commands:

```bash
gfortran -O6 -ffree-line-length-0 endpoint.f90 dopri5_no_rpar.f90 dopri5copy_no_rpar.f90 -o endpoint
```
or

```bash
ifort -xHost -O3 -w -qopt-report -prof-use -ipo endpoint.f90 dopri5_no_rpar.f90 dopri5copy_no_rpar.f90  -o endpoint
```

Then we need to just run perl file:

```bash
./plotendpoints.pl
```

or in order to save data we have:

```bash
./plotendpoints.pl > saved.txt
```

In case the libraries are not sourced:

```bash
which perl 
```
add the directory to cshrc:

```bash
open ~/.cshrc 
```

```bash
source ~/.cshrc
```

Useful commands:
If you suspended the jobs by (ctr+z instead of ctr+c), you can see them using:

```bash
fg
```
Then we can cancel them using 'ctr+c'



