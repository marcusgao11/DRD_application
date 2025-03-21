Steps to run:

1. install intel oneapi for fortran 

2. cd to lib directory and run make to build dependent lib
   fftw3.3.10, libtestlib

3. cd to project folder and run make

4. run the peoject by:  ./mcl6

5. view results in MATLAB:  configure the data_path and run: test_2.m 



Content:

1. Precipitation_dynamics: FIG. 14. 



Code structure:

Precipitation_dynamics:;
PhSolve: interface evolve of precipitation
diffusion_eqn: solute diffusion in the bulk
VeloSolve_combined: velocity of NS equations
PreSolve_2: pressure solver
