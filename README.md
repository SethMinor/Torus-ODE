# Torus-ODE
> Repo for MATLAB material somehow related to the motion of vortices on curved surfaces. Some related resources:

1. [Repo link](https://github.com/SethMinor/Torus-ODE)

2. [Research advisor](https://carretero.sdsu.edu/publications/index.html#papers)

## *vortex_simulation.m*
> Driver code for vortex simulation; choose how many vortices you'd like to see the dynamics here, as well as there starting locations (isothermal coordinates).

## *jacobitheta1.m* and *jacobitheta1.m*
> Crude MATLAB function files to compute the first several terms (double precision is usually reached before a truncation of 15 or so terms) of the first Jacobi theta function and its derivative. 

## *vortex_velocity.m* and *complexpotential.m*
> MATLAB functions that computes the physical velocity of a set of vortices and their corresponding complex potential function, respectively, using a relation involving Jacobi theta functions is used (see [_Fetter et. al._](https://journals.aps.org/pra/pdf/10.1103/PhysRevA.101.053606?casa_token=Y-7DK7Ny6GYAAAAA%3A6d0WPKGSS2jhegscwXxLSe6u0O6XRoSd-A1o1ET2RzNMRYmkRlpXAkEkiH7Ydck_I-JDhGq016_pfQ) for more details). 
