# no-slip at four boundaries & \partial_n p=0 
# \partial_n\phi = 0.
THINOBJ=mcl6.o inout.o finit.o dist.o interpolate_coef.o setup.o  global.o dissipation.o PhSolve_ori.o VeloSolve_combined.o\
VeloSolve_2_combined.o DCT.o ESOLP.o PreSolve_2.o PreSolve_init.o dct_fftw.o ch_bulk_energy.o\
CFT.o bcd.o matrix_tool.o particles_params.o particles.o particles_move.o reinint_by_PhSolve.o ph_mu_solve.o ch_bdr_condition.o\


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
#LIBSTATIC=~/testlib.a 
LIBSTATICPATH=-L$(LIBHOME)


LIBSTATIC=-ltestlib
OPT =  -r8 -O2 #-g # -qopenmp #-warn all #-gen-interfaces -warn interfaces# -openmp# -O  debug
#OPT =  -r8 -O2 -qopenmp # -gen-interfaces -warn interfaces #-stubs# -gen-interfaces -warn interfaces# -openmp# -O
FF77 =  ifort


#LIBSTATIC=-ltestlib_gfortran
#OPT = -fdefault-real-8 -g -Wall ##-DFULL_IMPLICIT
#FF77 = gfortran


all:  mcl6 
mcl6: $(THINOBJ)
	$(FF77) $(OPT) -o mcl6 $(THINOBJ) -L$(FFTW_LIBPATH)/lib $(FFTW_LIBOPT) $(LIBSTATICPATH) $(LIBSTATIC) $(LIBPATH1) $(OPT_Link1) 
mcl6.o: mcl6.f90 global.o particles_move.o 
	$(FF77)  $(OPT) -c mcl6.f90 
inout.o: inout.f90 
	$(FF77)  $(OPT) -c inout.f90
global.o: global.f90
	$(FF77)  $(OPT) -c global.f90

particles_params.o: particles_params.f90
	$(FF77)  $(OPT) -c particles_params.f90
	
particles.o: particles.f90 particles_params.o
	$(FF77)  $(OPT) -c particles.f90
particles_move.o: particles_move.f90 particles.o
	$(FF77)  $(OPT) -c particles_move.f90

finit.o: finit.f90 global.f90 particles.o
	$(FF77)  $(OPT) -c finit.f90
reinint_by_PhSolve.o: reinint_by_PhSolve.f90 global.f90
	$(FF77)  $(OPT) -c reinint_by_PhSolve.f90
dist.o: dist.f90 global.f90 particles_move.o
	$(FF77)  $(OPT) -c dist.f90
interpolate_coef.o: interpolate_coef.f90 global.f90
	$(FF77)  $(OPT) -c interpolate_coef.f90
setup.o: setup.f90 global.f90
	$(FF77)  $(OPT) -c setup.f90

PhSolve_ori.o: PhSolve_ori.f90 global.f90
	$(FF77)  $(OPT) -c PhSolve_ori.f90
ph_mu_solve.o: ph_mu_solve.f90 global.f90
	$(FF77)  $(OPT) -c ph_mu_solve.f90

matrix_tool.o: matrix_tool.f90 global.f90
	$(FF77)  $(OPT) -c matrix_tool.f90
PreSolve_init.o: PreSolve_init.f90
	$(FF77)  $(OPT) -I$(FFTW_LIBPATH)/include -c PreSolve_init.f90
PreSolve_2.o: PreSolve_2.f90 global.f90
	$(FF77)  $(OPT) -c PreSolve_2.f90
dct_fftw.o: dct_fftw.f90
	$(FF77)  $(OPT) -c dct_fftw.f90
bcd.o: bcd.f90 global.f90 
	$(FF77)  $(OPT) -c bcd.f90

VeloSolve_combined.o: VeloSolve_combined.f90 global.f90
	$(FF77)  $(OPT) -c VeloSolve_combined.f90
VeloSolve_2_combined.o: VeloSolve_2_combined.f90 global.f90
	$(FF77)  $(OPT) -c VeloSolve_2_combined.f90
dissipation.o: dissipation.f90 global.f90
	$(FF77)  $(OPT) -c dissipation.f90

ESOLP.o: ESOLP.f90
	$(FF77)  $(OPT) -c ESOLP.f90
DCT.o: DCT.f90
	$(FF77)  $(OPT) -c DCT.f90
CFT.o: CFT.f90
	$(FF77)  $(OPT) -c CFT.f90
ch_bulk_energy.o: ch_bulk_energy.f90
	$(FF77)  $(OPT) -c ch_bulk_energy.f90
ch_bdr_condition.o: ch_bdr_condition.f90
	$(FF77)  $(OPT) -c ch_bdr_condition.f90

clean:
	rm *.o *.mod mcl6 data data* srf fort.*
cleandata:
	rm data* srf fort.* 
cleano:
	rm *.o mcl6
