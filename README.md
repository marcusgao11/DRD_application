
# DRD-Micro （Diffuse-Resistance-Domain Approach）：A Unified Simulation Method of Microscale Multiphase Flows

## Steps to run:

1. install intel oneapi for fortran 

2. cd to lib directory and run make to build dependent lib
   fftw3.3.10, libtestlib

3. cd to project folder and run make

4. run the peoject by:  ./mcl6

5. view results in MATLAB:  configure the data_path and run: test_2.m 



Content: (Applications in the paper: [1])

0. two_phase_flow_complex_domain: FIG. 10.
1. Precipitation_dynamics: FIG. 14. 
2. active_matter_hydrodynamics: FIG. 12.
3. microfluidics_complex_domain: FIG. 11.


Code structure:

two_phase_flow_complex_domain:
ph_mu_solve: interface evolve of two-phase flows
VeloSolve_combined: velocity of NS equations
PreSolve_2: pressure solver
particles_move: 
particles_params: 
particles: solvers for particles and fiber
finit: complex domain is specified here


Precipitation_dynamics:;
PhSolve: interface evolve of precipitation
diffusion_eqn: solute diffusion in the bulk
VeloSolve_combined: velocity of NS equations
PreSolve_2: pressure solver

active_matter_hydrodynamics
VeloSolve_combined: velocity of NS equations
PreSolve_2: pressure solver
particles_move: 
particles_params: 
particles: solvers for particles and fiber

microfluidics_complex_domain:
PhSolve: interface evolve of two-phase flows
VeloSolve_combined: velocity of NS equations
PreSolve_2: pressure solver
particles_move: 
particles_params: 
particles: solvers for particles and fiber
finit: read snake structure from bin file


## Example: two-phase flows with fibers in complex channel
To run the example:

    > cd DRD_application/microfluidics_complex_domain  
    > make  
    > ./mcl6  

Running times on single CPU intel i7-14700:  around 4 hours

Plot the results:  

    > python .\post_process\plot_movie.py  
    > python .\post_process\plot_snapshot.py  

Results are also presented here:
* [two-phase flows with fibers in complex channel](doc/two-phase-flows-channel.md);

## Reference
[1] Gao, M., Li, Z., & Xu, X. (2024). Simulation Method of Microscale Fluid-Structure Interactions: Diffuse-Resistance-Domain Approach. arXiv preprint arXiv:2407.14291.
