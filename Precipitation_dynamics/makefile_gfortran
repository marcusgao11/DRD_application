# no-slip at four boundaries & \partial_n p=0 
# \partial_n\phi = 0.
THINOBJ=mcl6.o inout.o finit.o dist.o interpolate_coef.o setup.o\
global.o dissipation.o PhSolve.o VeloSolve_combined.o\
DCT.o ESOLP.o PreSolve_2.o PreSolve_init.o\
CFT.o bcd.o matrix_tool.o diffusion_eqn.o


# home directory
LIBHOME=../lib

##machine dependent
LIBPATH=/opt/intel/oneapi/mkl/2022.0.2
LIBPATH1=-L$(LIBPATH)/lib/intel64 -I$(LIBPATH)/include
OPT_Link1= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread


# FFTW related
#FFTW_LIBPATH = $(LIBHOME)/libFFT
FFTW_LIBPATH = $(LIBHOME)/fftw3.3.10
FFTW_LIBOPT = -lfftw3

# home directory
LIBSTATICPATH=-L$(LIBHOME)

LIBSTATIC=-ltestlib

OPT =  -r8 -O2 #-openmp #-stubs# -gen-interfaces -warn interfaces# -openmp# -O
FF77 =  ifort


## debug mode
#OPT = -freal-4-real-8 -g -Wall ##-DFULL_IMPLICIT
#FF77 = gfortran

all:  mcl6 
mcl6: $(THINOBJ)
	$(FF77) $(OPT) -o mcl6 $(THINOBJ) -L$(FFTW_LIBPATH)/lib $(FFTW_LIBOPT) $(LIBSTATICPATH) $(LIBSTATIC) $(LIBPATH1) $(OPT_Link1) 
mcl6.o: mcl6.f90 global.o
	$(FF77)  $(OPT) -c mcl6.f90 
inout.o: inout.f90 
	$(FF77)  $(OPT) -c inout.f90
global.o: global.f90
	$(FF77)  $(OPT) -c global.f90
finit.o: finit.f90 global.f90
	$(FF77)  $(OPT) -c finit.f90
dist.o: dist.f90 global.f90
	$(FF77)  $(OPT) -c dist.f90
interpolate_coef.o: interpolate_coef.f90 global.f90
	$(FF77)  $(OPT) -c interpolate_coef.f90
setup.o: setup.f90 global.f90
	$(FF77)  $(OPT) -c setup.f90
PhSolve.o: PhSolve.f90 global.f90
	$(FF77)  $(OPT) -c PhSolve.f90

diffusion_eqn.o: diffusion_eqn.f90 global.f90
	$(FF77)  $(OPT) -c diffusion_eqn.f90
matrix_tool.o: matrix_tool.f90 global.f90
	$(FF77)  $(OPT) -c matrix_tool.f90
PreSolve_init.o: PreSolve_init.f90
	$(FF77)  $(OPT) -I$(FFTW_LIBPATH)/include -c PreSolve_init.f90
PreSolve_2.o: PreSolve_2.f90 global.f90
	$(FF77)  $(OPT) -c PreSolve_2.f90
bcd.o: bcd.f90 global.f90 
	$(FF77)  $(OPT) -c bcd.f90
VeloSolve_combined.o: VeloSolve_combined.f90 global.f90
	$(FF77)  $(OPT) -c VeloSolve_combined.f90
dissipation.o: dissipation.f90 global.f90
	$(FF77)  $(OPT) -c dissipation.f90
ESOLP.o: ESOLP.f90
	$(FF77)  $(OPT) -c ESOLP.f90
DCT.o: DCT.f90
	$(FF77)  $(OPT) -c DCT.f90
CFT.o: CFT.f90
	$(FF77)  $(OPT) -c CFT.f90

clean:
	rm *.o *.mod mcl6 data data* srf fort.*
cleandata:
	rm data* srf fort.* 
cleano:
	rm *.o mcl6
